"""Common contract for every structural resource.

The three resources BioSeqInsight talks to are completely different: ESM Atlas
runs live inference on a POSTed sequence, AlphaFold DB serves a precomputed
file keyed by UniProt accession, and RCSB serves experimental coordinates
keyed by PDB identifier or found by sequence search. v1.0 wired each of them
straight into the GUI with its own ad-hoc error handling and return format.

v2.0 puts them behind one interface that returns one result type. That is
what makes it possible to compare the resources scientifically, to fall back
from one to another, and to report all three in a single exported table.
"""

from __future__ import annotations

import re
from abc import ABC, abstractmethod
from datetime import datetime, timezone
from pathlib import Path

from .. import __version__
from ..config.settings import Settings
from ..models.results import RetrievalStatus, StructureResult
from ..services.cache import ResponseCache, make_key
from ..services.http_client import HttpClient
from ..services.logging_setup import get_logger, sequence_reference

UNIPROT_PATTERN = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$",
    re.IGNORECASE,
)
PDB_ID_PATTERN = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def is_uniprot_accession(text: str) -> bool:
    return bool(UNIPROT_PATTERN.match(normalise_accession(text)))


def is_pdb_id(text: str) -> bool:
    token = (text or "").strip()
    return bool(PDB_ID_PATTERN.match(token))


def normalise_accession(text: str) -> str:
    """Extract a bare UniProt accession from the forms users actually paste.

    Handles ``P24941``, the FASTA header form ``sp|P24941|CDK2_HUMAN`` (where
    the accession is the middle field, not the last), and AlphaFold model
    names such as ``AF-P69905-F1``. Falls back to the raw token when nothing
    accession-shaped is present, so the caller can report it back to the user.
    """
    raw = (text or "").strip()
    if not raw:
        return ""
    token = raw.split()[0].upper()
    candidates = [part for part in token.replace("-", "|").split("|") if part]
    for candidate in candidates:
        if UNIPROT_PATTERN.match(candidate):
            return candidate
    return token


def classify_query(text: str) -> str:
    """Return ``"uniprot"``, ``"pdb_id"``, ``"sequence"`` or ``"unknown"``."""
    raw = (text or "").strip()
    if not raw:
        return "unknown"
    token = raw.split()[0]
    if is_uniprot_accession(token) and len(raw.split()) == 1 and len(raw) <= 15:
        return "uniprot"
    if is_pdb_id(token) and len(raw) <= 6:
        return "pdb_id"
    if len(raw) >= 10:
        return "sequence"
    return "unknown"


class StructureProvider(ABC):
    """Base class for ESM Atlas, AlphaFold DB and RCSB PDB providers."""

    #: short identifier used in results, cache namespaces and benchmark tables
    name: str = "provider"
    #: whether this provider returns predicted (as opposed to experimental) models
    predicted: bool = True

    def __init__(
        self,
        client: HttpClient,
        settings: Settings,
        cache: ResponseCache | None = None,
    ) -> None:
        self.client = client
        self.settings = settings
        self.cache = cache
        self.log = get_logger(self.name)

    # -- capabilities -------------------------------------------------------

    @abstractmethod
    def supports(self, query_type: str) -> bool:
        """Can this provider answer a query of this kind?"""

    @abstractmethod
    def retrieve(
        self,
        query: str,
        *,
        query_type: str | None = None,
        sequence: str | None = None,
        accession: str | None = None,
        refresh: bool = False,
    ) -> StructureResult:
        """Fetch a structure and return a populated :class:`StructureResult`."""

    # -- shared helpers -----------------------------------------------------

    def _new_result(self, query: str, query_type: str) -> StructureResult:
        return StructureResult(
            source=self.name,
            query=query,
            query_type=query_type,
            status=RetrievalStatus.SKIPPED,
            experimental=not self.predicted,
            retrieved_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
            software_version=__version__,
        )

    def _cache_key(self, query: str, **params) -> str:
        return make_key(self.name, query, **params)

    def _cache_get(self, key: str, refresh: bool) -> str | None:
        if self.cache is None:
            return None
        hit = self.cache.get(self.name, key, refresh=refresh)
        return hit[0] if hit else None

    def _cache_put(self, key: str, payload: str, query: str, url: str | None) -> None:
        if self.cache is None or not self.cache.enabled:
            return
        try:
            self.cache.put(self.name, key, payload, query=query, source=self.name, url=url)
        except Exception as exc:  # pragma: no cover - cache must never break retrieval
            self.log.warning("Could not cache %s response: %s", self.name, exc)

    def _save(self, pdb_text: str, filename: str) -> str:
        directory = Path(self.settings.download_dir)
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / filename
        path.write_text(pdb_text, encoding="utf-8")
        return str(path)

    def _log_query(self, operation: str, query: str, query_type: str) -> None:
        reference = sequence_reference(query) if query_type == "sequence" else query
        self.log.info("%s %s (%s)", operation, reference, query_type)


__all__ = [
    "PDB_ID_PATTERN",
    "StructureProvider",
    "UNIPROT_PATTERN",
    "classify_query",
    "is_pdb_id",
    "is_uniprot_accession",
    "normalise_accession",
]
