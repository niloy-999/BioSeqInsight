"""The structure subsystem's front door.

:class:`StructureManager` accepts whatever the user has (a sequence, a UniProt
accession or a PDB identifier), decides which resources can answer, queries
them in a configurable order, validates every returned structure against the
query, and returns a ranked list of :class:`StructureResult` objects.

The fallback rule that matters
------------------------------
Falling back is only useful if it never quietly substitutes a different
molecule. The manager therefore stops early **only** on an M4 (exact)
result. Anything less than exact is kept, the remaining resources are still
tried, and the user is shown every candidate with its mapping level. A
warning is attached whenever the best available answer falls short of exact
identity, which is precisely the situation the v1.0 benchmark hit 7 times out
of 25 and reported without comment.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

from ..config.settings import Settings, get_settings
from ..core.alphabet import clean_protein
from ..models.results import MappingLevel, RetrievalStatus, StructureResult
from ..services.cache import NullCache, ResponseCache
from ..services.errors import BioSeqInsightError
from ..services.http_client import HttpClient, RetryPolicy
from ..services.logging_setup import get_logger
from .alphafold import AlphaFoldProvider
from .base import classify_query, normalise_accession
from .esmatlas import ESMAtlasProvider
from .rcsb import RCSBProvider
from .uniprot import UniProtClient
from .validation import no_result_report, validate_structure

logger = get_logger("structures")


@dataclass
class StructureQuery:
    """A resolved query: what the user typed plus everything inferred from it."""

    raw: str
    query_type: str
    sequence: str | None = None
    accession: str | None = None
    pdb_id: str | None = None
    notes: list[str] = field(default_factory=list)


@dataclass
class StructureSearchOutcome:
    """Everything the manager learned about one query."""

    query: StructureQuery
    results: list[StructureResult] = field(default_factory=list)
    elapsed_s: float = 0.0
    warnings: list[str] = field(default_factory=list)

    @property
    def best(self) -> StructureResult | None:
        usable = [r for r in self.results if r.ok]
        if not usable:
            return None
        return max(
            usable,
            key=lambda r: (
                r.mapping_level.rank,
                r.validation.identity_percent or 0.0 if r.validation else 0.0,
                1 if r.experimental else 0,
            ),
        )

    @property
    def best_mapping(self) -> MappingLevel:
        best = self.best
        return best.mapping_level if best else MappingLevel.M0_NO_RESULT

    def to_dict(self) -> dict:
        return {
            "query": self.query.raw,
            "query_type": self.query.query_type,
            "accession": self.query.accession,
            "elapsed_s": self.elapsed_s,
            "best_mapping": self.best_mapping.value,
            "best_mapping_label": self.best_mapping.label,
            "warnings": list(self.warnings),
            "results": [r.to_dict() for r in self.results],
        }


class StructureManager:
    """Coordinates every structural resource behind one unified interface."""

    def __init__(
        self,
        settings: Settings | None = None,
        client: HttpClient | None = None,
        cache: ResponseCache | None = None,
    ) -> None:
        self.settings = settings or get_settings()
        self.cache = cache if cache is not None else self._default_cache(self.settings)
        self.client = client or HttpClient(
            policy=RetryPolicy(
                max_attempts=self.settings.max_attempts,
                backoff_base_s=self.settings.backoff_base_s,
                backoff_max_s=self.settings.backoff_max_s,
                jitter=self.settings.jitter,
            ),
            timeout_s=self.settings.http_timeout_s,
        )
        self.providers = {
            "rcsb": RCSBProvider(self.client, self.settings, self.cache),
            "alphafold": AlphaFoldProvider(self.client, self.settings, self.cache),
            "esmatlas": ESMAtlasProvider(self.client, self.settings, self.cache),
        }
        self.uniprot = UniProtClient(self.client, self.settings, self.cache)

    @staticmethod
    def _default_cache(settings: Settings) -> ResponseCache:
        if not settings.cache_enabled:
            return NullCache()
        return ResponseCache(settings.cache_dir, ttl_s=settings.cache_ttl_s)

    # -- query resolution ---------------------------------------------------

    def resolve(self, text: str, refresh: bool = False) -> StructureQuery:
        """Work out what the user gave us and fill in what can be derived.

        A UniProt accession is expanded to its canonical sequence when the
        network allows, because without the sequence no retrieved structure
        can be validated and every result would be capped at M1.
        """
        raw = (text or "").strip()
        kind = classify_query(raw)
        query = StructureQuery(raw=raw, query_type=kind)

        if kind == "uniprot":
            query.accession = normalise_accession(raw)
            if not self.settings.offline:
                try:
                    query.sequence = self.uniprot.fetch_sequence(query.accession, refresh=refresh)
                except BioSeqInsightError as exc:
                    query.notes.append(
                        f"Could not fetch the canonical sequence for {query.accession} "
                        f"({exc}). Retrieved structures cannot be identity-checked, so "
                        "results will be reported as unverified (M1)."
                    )
        elif kind == "pdb_id":
            query.pdb_id = raw.upper()
        elif kind == "sequence":
            query.sequence = clean_protein(raw)
        return query

    # -- retrieval ----------------------------------------------------------

    def retrieve(
        self,
        text: str,
        *,
        providers: list[str] | None = None,
        stop_on_exact: bool = True,
        refresh: bool = False,
        validate: bool | None = None,
    ) -> StructureSearchOutcome:
        """Query the configured resources and validate every hit."""
        started = time.perf_counter()
        query = self.resolve(text, refresh=refresh)
        outcome = StructureSearchOutcome(query=query, warnings=list(query.notes))
        should_validate = self.settings.validate_structures if validate is None else validate

        if query.query_type == "unknown":
            failed = StructureResult(
                source="manager",
                query=query.raw,
                query_type="unknown",
                status=RetrievalStatus.INVALID_INPUT,
                error=(
                    "Input was not recognised as a protein sequence, a UniProt accession "
                    "or a four-character PDB identifier."
                ),
                validation=no_result_report("Input could not be interpreted."),
            )
            outcome.results.append(failed)
            outcome.elapsed_s = time.perf_counter() - started
            return outcome

        order = providers or list(self.settings.fallback_order)
        for name in order:
            provider = self.providers.get(name)
            if provider is None:
                outcome.warnings.append(f"Unknown provider {name!r} was skipped.")
                continue
            effective_type = self._effective_type(name, query)
            if effective_type is None:
                continue

            try:
                result = provider.retrieve(
                    query.raw,
                    query_type=effective_type,
                    sequence=query.sequence,
                    accession=query.accession,
                    refresh=refresh,
                )
            except Exception as exc:  # pragma: no cover - providers convert their own errors
                logger.exception("Provider %s raised unexpectedly", name)
                result = StructureResult(
                    source=name,
                    query=query.raw,
                    query_type=effective_type,
                    status=RetrievalStatus.SERVICE_ERROR,
                    error=f"{type(exc).__name__}: {exc}",
                )

            if result.ok and should_validate:
                result.validation = validate_structure(
                    self._read_back(result),
                    query_sequence=query.sequence,
                    query_accession=query.accession,
                    settings=self.settings,
                    extra_accessions=result.cross_references,
                )
            elif not result.ok and result.validation is None:
                result.validation = no_result_report(result.error or "")

            outcome.results.append(result)
            if stop_on_exact and result.ok and result.mapping_level is MappingLevel.M4_EXACT:
                break

        outcome.elapsed_s = time.perf_counter() - started
        self._attach_warnings(outcome)
        return outcome

    def _effective_type(self, provider: str, query: StructureQuery) -> str | None:
        """Decide how a provider should be queried for this input, or skip it."""
        if provider == "alphafold":
            return "uniprot" if query.accession else None
        if provider == "esmatlas":
            return "sequence" if query.sequence else None
        if provider == "rcsb":
            if query.pdb_id:
                return "pdb_id"
            if query.sequence:
                return "sequence"
            return None
        return None  # pragma: no cover - guarded by settings.validate()

    def _read_back(self, result: StructureResult) -> str:
        """Read the saved coordinate file for validation."""
        if not result.download_path:
            return ""
        try:
            with open(result.download_path, encoding="utf-8", errors="replace") as handle:
                return handle.read()
        except OSError as exc:  # pragma: no cover - filesystem failure
            logger.warning("Cannot re-read %s for validation: %s", result.download_path, exc)
            return ""

    def _attach_warnings(self, outcome: StructureSearchOutcome) -> None:
        best = outcome.best
        if best is None:
            outcome.warnings.append(
                "No resource returned a usable structure for this query. "
                "The requested molecule may have no deposited or predicted model, "
                "or the services may be unavailable."
            )
            return
        level = best.mapping_level
        if level is MappingLevel.M4_EXACT:
            return
        if level is MappingLevel.M3_HIGH_CONFIDENCE:
            outcome.warnings.append(
                "The best available structure is a high-confidence sequence match but its "
                "biological identity was not confirmed against the requested accession."
            )
        elif level is MappingLevel.M2_SEQUENCE_MATCH:
            outcome.warnings.append(
                "The best available structure is a related record, not the requested "
                "molecule. Do not report it as the structure of your protein."
            )
        else:
            outcome.warnings.append(
                "A structure was retrieved but could not be identity-checked. "
                "Supply a sequence or accession to enable validation."
            )

    # -- direct helpers used by the GUI -------------------------------------

    def fetch_pdb_id(self, pdb_id: str, refresh: bool = False) -> StructureResult:
        return self.providers["rcsb"].fetch_entry(pdb_id, refresh=refresh)

    def fetch_alphafold(self, accession: str, refresh: bool = False) -> StructureResult:
        return self.providers["alphafold"].retrieve(
            accession,
            query_type="uniprot",
            accession=normalise_accession(accession),
            refresh=refresh,
        )

    def fold_sequence(self, sequence: str, refresh: bool = False) -> StructureResult:
        return self.providers["esmatlas"].retrieve(
            sequence, query_type="sequence", sequence=sequence, refresh=refresh
        )

    def cache_stats(self) -> dict:
        return self.cache.stats()


__all__ = ["StructureManager", "StructureQuery", "StructureSearchOutcome"]
