"""ESM Atlas provider: live remote ESMFold inference.

BioSeqInsight does not run ESMFold. It POSTs a sequence to the public ESM
Atlas fold endpoint and stores whatever coordinates come back. That service
is a shared demo endpoint and fails often: in the v1.0 benchmark it returned
HTTP 504 for 9 of 25 requests. This provider therefore does three things v1.0
did not:

* it goes through the retrying :class:`~bioseqinsight.services.http_client.HttpClient`,
  so transient 5xx responses get a bounded number of further chances and each
  attempt is recorded;
* it validates that a 200 response actually contains coordinates, because a
  gateway error page returned with status 200 is not a structure;
* it refuses sequences above a configurable length rather than sending a
  request that is known to time out, and says so instead of failing silently.
"""

from __future__ import annotations

import time
from datetime import datetime, timezone

from ..config.settings import ESMATLAS_FOLD_URL
from ..core.alphabet import SequenceError, validate_protein
from ..models.results import RetrievalStatus, StructureResult
from .base import StructureProvider
from .pdbio import looks_like_pdb, mean_confidence, parse_pdb


class ESMAtlasProvider(StructureProvider):
    """Remote ESMFold prediction via the ESM Atlas fold API."""

    name = "esmatlas"
    predicted = True

    def supports(self, query_type: str) -> bool:
        return query_type == "sequence"

    def retrieve(
        self,
        query: str,
        *,
        query_type: str | None = None,
        sequence: str | None = None,
        accession: str | None = None,
        refresh: bool = False,
    ) -> StructureResult:
        result = self._new_result(query, query_type or "sequence")
        raw = sequence or query

        try:
            protein = validate_protein(raw)
        except SequenceError as exc:
            result.status = RetrievalStatus.INVALID_INPUT
            result.error = str(exc)
            return result

        if len(protein) > self.settings.esm_max_length:
            result.status = RetrievalStatus.INVALID_INPUT
            result.error = (
                f"Sequence is {len(protein)} aa. The public ESM Atlas endpoint is "
                f"configured here for at most {self.settings.esm_max_length} aa because "
                "longer requests reliably time out. Use AlphaFold DB or RCSB for this "
                "protein, or raise esm_max_length in the settings."
            )
            return result

        if self.settings.offline:
            result.status = RetrievalStatus.SKIPPED
            result.error = "Offline mode: live ESM Atlas folding was not attempted."
            return result

        key = self._cache_key(protein, endpoint=ESMATLAS_FOLD_URL)
        started = time.perf_counter()
        cached = self._cache_get(key, refresh)
        if cached is not None:
            self._finalise(result, cached, protein, from_cache=True)
            result.response_time_s = time.perf_counter() - started
            return result

        self._log_query("ESM Atlas fold", protein, "sequence")
        outcome = self.client.post(
            ESMATLAS_FOLD_URL,
            data=protein.encode("utf-8"),
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            timeout_s=self.settings.esm_timeout_s,
            expect_text="ATOM",
        )
        result.attempts = outcome.attempts
        result.response_time_s = outcome.total_elapsed_s

        if outcome.ok and outcome.response is not None:
            self._cache_put(key, outcome.response.text, protein, ESMATLAS_FOLD_URL)
            self._finalise(result, outcome.response.text, protein, from_cache=False)
            return result

        last = outcome.attempts[-1] if outcome.attempts else None
        timed_out = bool(last and last.error and last.error.startswith("timeout"))
        result.status = RetrievalStatus.TIMEOUT if timed_out else RetrievalStatus.SERVICE_ERROR
        result.structure_url = ESMATLAS_FOLD_URL
        result.error = (
            f"ESM Atlas did not return coordinates after {len(outcome.attempts)} attempt(s). "
            f"Last error: {last.error if last else 'unknown'}. "
            "A 504 here means the public service timed out, not that the sequence is invalid."
        )
        return result

    def _finalise(
        self, result: StructureResult, pdb_text: str, protein: str, from_cache: bool
    ) -> None:
        if not looks_like_pdb(pdb_text):
            result.status = RetrievalStatus.SERVICE_ERROR
            result.error = "ESM Atlas returned a payload without ATOM records."
            return
        parsed = parse_pdb(pdb_text)
        confidence, kind = mean_confidence(pdb_text, is_predicted=True)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
        filename = f"esmfold_{len(protein)}aa_{stamp}.pdb"
        result.status = RetrievalStatus.OK
        result.identifier = f"ESMFold:{len(protein)}aa"
        result.structure_url = ESMATLAS_FOLD_URL
        result.download_path = self._save(pdb_text, filename)
        result.sequence = parsed.longest_chain.sequence if parsed.longest_chain else None
        result.confidence = confidence
        result.confidence_kind = kind
        result.from_cache = from_cache


__all__ = ["ESMAtlasProvider"]
