"""UniProt lookups.

UniProt is not a structure provider, but it is what makes validation possible:
to decide whether a retrieved PDB chain really is the protein behind
accession P24941, BioSeqInsight needs P24941's canonical sequence. This
client fetches that sequence (and, when asked, the entry's cross-referenced
PDB identifiers) through the public REST API, with the same retry policy and
cache as every other external call.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

from ..config.settings import (
    UNIPROT_FASTA_URL,
    UNIPROT_JSON_URL,
    Settings,
)
from ..core.alphabet import clean_protein, parse_fasta
from ..services.cache import ResponseCache, make_key
from ..services.errors import NotFoundError, ServiceUnavailable
from ..services.http_client import HttpClient
from ..services.logging_setup import get_logger
from .base import normalise_accession

logger = get_logger("uniprot")


@dataclass
class UniProtEntry:
    """The subset of a UniProt record BioSeqInsight uses."""

    accession: str
    sequence: str
    name: str = ""
    organism: str = ""
    length: int = 0
    pdb_ids: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "accession": self.accession,
            "name": self.name,
            "organism": self.organism,
            "length": self.length,
            "pdb_ids": list(self.pdb_ids),
        }


class UniProtClient:
    """Fetch canonical sequences and cross-references from UniProtKB."""

    name = "uniprot"

    def __init__(
        self,
        client: HttpClient,
        settings: Settings,
        cache: ResponseCache | None = None,
    ) -> None:
        self.client = client
        self.settings = settings
        self.cache = cache

    def fetch_sequence(self, accession: str, refresh: bool = False) -> str:
        """Return the canonical sequence for an accession."""
        acc = normalise_accession(accession)
        if not acc:
            raise NotFoundError("Empty UniProt accession.")
        key = make_key(self.name, acc, kind="fasta")
        cached = self._cache_get(key, refresh)
        if cached is None:
            url = UNIPROT_FASTA_URL.format(accession=acc)
            outcome = self.client.get(url, expect_text=">")
            if not outcome.ok or outcome.response is None:
                if outcome.status_code == 404:
                    raise NotFoundError(f"UniProt has no entry for {acc}.")
                raise ServiceUnavailable("UniProt", len(outcome.attempts), str(outcome.error or ""))
            cached = outcome.response.text
            self._cache_put(key, cached, acc, url)
        records = parse_fasta(cached)
        if not records:
            raise NotFoundError(f"UniProt returned no sequence for {acc}.")
        return clean_protein(records[0].sequence)

    def fetch_entry(self, accession: str, refresh: bool = False) -> UniProtEntry:
        """Return sequence plus name, organism and cross-referenced PDB ids."""
        acc = normalise_accession(accession)
        key = make_key(self.name, acc, kind="json")
        payload = self._cache_get(key, refresh)
        if payload is None:
            url = UNIPROT_JSON_URL.format(accession=acc)
            outcome = self.client.get(url, expect_text="primaryAccession")
            if not outcome.ok or outcome.response is None:
                if outcome.status_code == 404:
                    raise NotFoundError(f"UniProt has no entry for {acc}.")
                raise ServiceUnavailable("UniProt", len(outcome.attempts), str(outcome.error or ""))
            payload = outcome.response.text
            self._cache_put(key, payload, acc, url)
        return self._parse_entry(acc, payload)

    @staticmethod
    def _parse_entry(accession: str, payload: str) -> UniProtEntry:
        try:
            data = json.loads(payload)
        except json.JSONDecodeError as exc:
            raise NotFoundError(f"Malformed UniProt payload for {accession}: {exc}") from exc
        sequence = clean_protein((data.get("sequence") or {}).get("value", ""))
        organism = (data.get("organism") or {}).get("scientificName", "")
        description = (
            ((data.get("proteinDescription") or {}).get("recommendedName") or {}).get("fullName")
            or {}
        ).get("value", "")
        pdb_ids = [
            ref.get("id", "")
            for ref in data.get("uniProtKBCrossReferences", [])
            if ref.get("database") == "PDB" and ref.get("id")
        ]
        return UniProtEntry(
            accession=data.get("primaryAccession", accession),
            sequence=sequence,
            name=description,
            organism=organism,
            length=len(sequence),
            pdb_ids=pdb_ids,
        )

    # -- cache helpers ------------------------------------------------------

    def _cache_get(self, key: str, refresh: bool) -> str | None:
        if self.cache is None:
            return None
        hit = self.cache.get(self.name, key, refresh=refresh)
        return hit[0] if hit else None

    def _cache_put(self, key: str, payload: str, query: str, url: str) -> None:
        if self.cache is None or not self.cache.enabled:
            return
        try:
            self.cache.put(self.name, key, payload, query=query, source=self.name, url=url)
        except Exception as exc:  # pragma: no cover - cache must never break retrieval
            logger.warning("Could not cache UniProt response: %s", exc)


__all__ = ["UniProtClient", "UniProtEntry"]
