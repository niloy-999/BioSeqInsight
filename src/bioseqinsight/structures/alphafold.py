"""AlphaFold Protein Structure Database provider.

AlphaFold DB serves precomputed models keyed by UniProt accession. This
provider queries the public prediction API first
(``/api/prediction/{accession}``) to discover the current model file URL,
rather than guessing a model-version number in the file path.

That distinction matters and was the subject of a real defect: an earlier
version of this provider hard-coded candidate versions (``v4``, then ``v3``)
directly into the file URL. EBI periodically republishes AlphaFold DB under a
new model version, and once the database moved past the hard-coded versions,
every single accession returned HTTP 404 from the file endpoint — which is
*indistinguishable from a genuine coverage gap* unless something checks the
authoritative source first. The prediction API returns the file URL that is
current right now, so this provider no longer needs to know or guess a
version number at all.

Two failure modes matter and are handled explicitly:

* **404 from the prediction API is a genuine coverage gap.** AlphaFold DB has
  no model for some entries (in the v1.0 benchmark, T4 lysozyme P00720
  returned 404). This must not be retried or reported as a service failure,
  because doing so inflates the apparent unreliability of the resource.
* **A malformed or empty API response, or a file download that fails after
  the API named it, is a service failure**, not a coverage gap, and is
  reported and retried as such.
"""

from __future__ import annotations

import json
import time

from ..config.settings import ALPHAFOLD_API_URL
from ..models.results import RetrievalStatus, StructureResult
from .base import StructureProvider, normalise_accession
from .pdbio import looks_like_pdb, mean_confidence, parse_pdb


class AlphaFoldProvider(StructureProvider):
    """Retrieve a precomputed AlphaFold DB model by UniProt accession."""

    name = "alphafold"
    predicted = True

    def supports(self, query_type: str) -> bool:
        return query_type == "uniprot"

    def retrieve(
        self,
        query: str,
        *,
        query_type: str | None = None,
        sequence: str | None = None,
        accession: str | None = None,
        refresh: bool = False,
    ) -> StructureResult:
        result = self._new_result(query, query_type or "uniprot")
        acc = normalise_accession(accession or query)

        if not acc:
            result.status = RetrievalStatus.INVALID_INPUT
            result.error = (
                "AlphaFold DB is keyed by UniProt accession and none was supplied. "
                "Example: P69905."
            )
            return result

        if self.settings.offline:
            result.status = RetrievalStatus.SKIPPED
            result.error = "Offline mode: AlphaFold DB was not queried."
            return result

        started = time.perf_counter()
        attempts_all: list = []
        api_url = ALPHAFOLD_API_URL.format(accession=acc)
        api_key = self._cache_key(acc, kind="api")
        cached_api = self._cache_get(api_key, refresh)

        if cached_api is not None:
            try:
                entries = json.loads(cached_api)
            except ValueError:
                entries = None
        else:
            self._log_query("AlphaFold DB prediction lookup", acc, "uniprot")
            outcome = self.client.get(api_url, expect_text="pdbUrl")
            attempts_all.extend(outcome.attempts)

            if outcome.status_code == 404:
                result.attempts = attempts_all
                result.response_time_s = time.perf_counter() - started
                result.structure_url = api_url
                result.status = RetrievalStatus.NOT_FOUND
                result.error = (
                    f"AlphaFold DB holds no model for {acc}. This is a coverage gap in "
                    "the database, not a service failure."
                )
                return result

            if not outcome.ok or outcome.response is None:
                result.attempts = attempts_all
                result.response_time_s = time.perf_counter() - started
                result.structure_url = api_url
                self._service_failure(result, attempts_all, "prediction lookup", acc)
                return result

            try:
                entries = json.loads(outcome.response.text)
            except ValueError:
                entries = None
            else:
                self._cache_put(api_key, outcome.response.text, acc, api_url)

        if not entries:
            result.attempts = attempts_all
            result.response_time_s = time.perf_counter() - started
            result.structure_url = api_url
            if entries is None:
                result.status = RetrievalStatus.SERVICE_ERROR
                result.error = "AlphaFold DB returned a response that was not valid JSON."
            else:
                result.status = RetrievalStatus.NOT_FOUND
                result.error = (
                    f"AlphaFold DB holds no model for {acc}. This is a coverage gap in "
                    "the database, not a service failure."
                )
            return result

        entry = entries[0]
        pdb_url = entry.get("pdbUrl")
        entry_id = entry.get("entryId") or f"AF-{acc}-F1"
        if not pdb_url:
            result.attempts = attempts_all
            result.response_time_s = time.perf_counter() - started
            result.status = RetrievalStatus.SERVICE_ERROR
            result.error = f"AlphaFold DB prediction entry for {acc} had no pdbUrl field."
            return result

        file_key = self._cache_key(acc, kind="file", entry=entry_id)
        cached_file = self._cache_get(file_key, refresh)
        if cached_file is not None:
            result.attempts = attempts_all
            self._finalise(result, cached_file, entry_id, pdb_url, from_cache=True)
            result.response_time_s = time.perf_counter() - started
            return result

        self._log_query("AlphaFold DB model download", entry_id, "uniprot")
        outcome = self.client.get(pdb_url, expect_text="ATOM")
        attempts_all.extend(outcome.attempts)
        result.attempts = attempts_all
        result.response_time_s = time.perf_counter() - started

        if not outcome.ok or outcome.response is None:
            result.structure_url = pdb_url
            self._service_failure(result, attempts_all, "model download", acc)
            return result

        self._cache_put(file_key, outcome.response.text, acc, pdb_url)
        self._finalise(result, outcome.response.text, entry_id, pdb_url, from_cache=False)
        return result

    def _service_failure(self, result, attempts_all, stage: str, accession: str) -> None:
        last = attempts_all[-1] if attempts_all else None
        timed_out = bool(last and last.error and last.error.startswith("timeout"))
        result.status = RetrievalStatus.TIMEOUT if timed_out else RetrievalStatus.SERVICE_ERROR
        result.error = (
            f"AlphaFold DB {stage} for {accession} failed after {len(attempts_all)} "
            f"attempt(s). Last error: {last.error if last else 'unknown'}"
        )

    def _finalise(
        self,
        result: StructureResult,
        pdb_text: str,
        entry_id: str,
        url: str,
        from_cache: bool,
    ) -> None:
        if not looks_like_pdb(pdb_text):
            result.status = RetrievalStatus.SERVICE_ERROR
            result.error = "AlphaFold DB returned a payload without ATOM records."
            return
        parsed = parse_pdb(pdb_text)
        confidence, kind = mean_confidence(pdb_text, is_predicted=True)
        result.status = RetrievalStatus.OK
        result.identifier = entry_id
        result.structure_url = url
        result.download_path = self._save(pdb_text, f"{entry_id}.pdb")
        result.sequence = parsed.longest_chain.sequence if parsed.longest_chain else None
        result.confidence = confidence
        result.confidence_kind = kind
        result.from_cache = from_cache


__all__ = ["AlphaFoldProvider"]
