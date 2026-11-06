"""RCSB PDB provider: experimental coordinates by identifier or by sequence.

Two access routes are implemented and, importantly, kept distinct in the
result:

``pdb_id``
    A direct, unambiguous download. Identity is whatever the entry says it is.

``sequence``
    A similarity search against the PDB sequence index. This is the route
    that produced the most misleading numbers in v1.0: the search returned a
    hit for every query, which looks like a 100% success rate until you check
    whether each hit is the right molecule. Here the search returns a ranked
    list, each candidate's polymer entities are fetched so their UniProt
    cross-references are known, and the caller validates identity before
    anything is presented as a match.
"""

from __future__ import annotations

import json
import time

from ..config.settings import (
    RCSB_DOWNLOAD_URL,
    RCSB_ENTRY_URL,
    RCSB_POLYMER_URL,
    RCSB_SEARCH_URL,
)
from ..core.alphabet import clean_protein
from ..models.results import RetrievalStatus, StructureResult
from .base import StructureProvider, is_pdb_id
from .pdbio import looks_like_pdb, mean_confidence, parse_pdb


class RCSBProvider(StructureProvider):
    """Experimental structures from the RCSB Protein Data Bank."""

    name = "rcsb"
    predicted = False

    def supports(self, query_type: str) -> bool:
        return query_type in ("pdb_id", "sequence", "uniprot")

    # -- public API ---------------------------------------------------------

    def retrieve(
        self,
        query: str,
        *,
        query_type: str | None = None,
        sequence: str | None = None,
        accession: str | None = None,
        refresh: bool = False,
    ) -> StructureResult:
        kind = query_type or ("pdb_id" if is_pdb_id(query) else "sequence")
        if self.settings.offline:
            result = self._new_result(query, kind)
            result.status = RetrievalStatus.SKIPPED
            result.error = "Offline mode: RCSB PDB was not queried."
            return result
        if kind == "pdb_id":
            return self.fetch_entry(query.strip().upper(), refresh=refresh)
        seq = clean_protein(sequence or query)
        return self.search_by_sequence(seq, accession=accession, refresh=refresh)

    def fetch_entry(
        self, pdb_id: str, *, refresh: bool = False, query_label: str | None = None
    ) -> StructureResult:
        """Download one experimental entry by its four-character identifier."""
        result = self._new_result(query_label or pdb_id, "pdb_id")
        pid = (pdb_id or "").strip().upper()
        if not is_pdb_id(pid):
            result.status = RetrievalStatus.INVALID_INPUT
            result.error = f"{pdb_id!r} is not a four-character PDB identifier (example: 1UBQ)."
            return result

        started = time.perf_counter()
        url = RCSB_DOWNLOAD_URL.format(pdb_id=pid)
        key = self._cache_key(pid, kind="coordinates")
        cached = self._cache_get(key, refresh)
        if cached is not None:
            self._finalise(result, cached, pid, url, from_cache=True)
            result.response_time_s = time.perf_counter() - started
            return result

        self._log_query("RCSB download", pid, "pdb_id")
        outcome = self.client.get(url, expect_text="ATOM")
        result.attempts = outcome.attempts
        result.response_time_s = outcome.total_elapsed_s
        if outcome.ok and outcome.response is not None:
            self._cache_put(key, outcome.response.text, pid, url)
            self._finalise(result, outcome.response.text, pid, url, from_cache=False)
            return result

        if outcome.status_code == 404:
            result.status = RetrievalStatus.NOT_FOUND
            result.error = f"RCSB has no entry {pid}."
            return result
        last = outcome.attempts[-1] if outcome.attempts else None
        timed_out = bool(last and last.error and last.error.startswith("timeout"))
        result.status = RetrievalStatus.TIMEOUT if timed_out else RetrievalStatus.SERVICE_ERROR
        result.error = f"RCSB download of {pid} failed: {last.error if last else 'unknown'}"
        return result

    def search_by_sequence(
        self,
        sequence: str,
        *,
        accession: str | None = None,
        refresh: bool = False,
    ) -> StructureResult:
        """Run a sequence search and download the best-scoring hit."""
        result = self._new_result(sequence, "sequence")
        seq = clean_protein(sequence)
        if len(seq) < 20:
            result.status = RetrievalStatus.INVALID_INPUT
            result.error = "RCSB sequence search needs at least 20 residues."
            return result

        started = time.perf_counter()
        try:
            hits = self.sequence_search_hits(seq, refresh=refresh)
        except RuntimeError as exc:
            result.status = RetrievalStatus.SERVICE_ERROR
            result.error = str(exc)
            result.response_time_s = time.perf_counter() - started
            return result

        if not hits:
            result.status = RetrievalStatus.NOT_FOUND
            result.error = (
                "RCSB sequence search returned no hit above the configured identity cutoff "
                f"({self.settings.rcsb_search_identity_cutoff:.0%})."
            )
            result.response_time_s = time.perf_counter() - started
            return result

        # Prefer a hit whose cross-references already contain the queried
        # accession. This alone converts many of the v1.0 "mapped to a
        # different protein" cases into correct retrievals, because the best
        # scoring hit is not always the one that is the requested molecule.
        ordered = hits
        if accession:
            wanted = accession.upper()
            ordered = sorted(
                hits,
                key=lambda h: (
                    0 if wanted in {a.upper() for a in h["accessions"]} else 1,
                    -h["score"],
                ),
            )

        chosen = ordered[0]
        entry = self.fetch_entry(chosen["pdb_id"], refresh=refresh, query_label=sequence)
        entry.query_type = "sequence"
        entry.response_time_s = time.perf_counter() - started
        if entry.ok:
            entry.identifier = chosen["pdb_id"]
            entry.structure_url = RCSB_DOWNLOAD_URL.format(pdb_id=chosen["pdb_id"])
            entry.cross_references = list(chosen["accessions"])
        return entry

    # -- search plumbing ----------------------------------------------------

    def sequence_search_hits(self, sequence: str, *, refresh: bool = False) -> list[dict]:
        """Return ranked search hits with their UniProt cross-references."""
        payload = json.dumps(self._search_payload(sequence)).encode("utf-8")
        key = self._cache_key(
            sequence,
            kind="search",
            cutoff=self.settings.rcsb_search_identity_cutoff,
            rows=self.settings.rcsb_search_max_hits,
        )
        raw = self._cache_get(key, refresh)
        if raw is None:
            self._log_query("RCSB sequence search", sequence, "sequence")
            outcome = self.client.post(
                RCSB_SEARCH_URL,
                data=payload,
                headers={"Content-Type": "application/json"},
            )
            if outcome.status_code == 204:
                return []
            if not outcome.ok or outcome.response is None:
                raise RuntimeError(
                    "RCSB sequence search failed after "
                    f"{len(outcome.attempts)} attempt(s): {outcome.error}"
                )
            raw = outcome.response.text
            self._cache_put(key, raw, sequence, RCSB_SEARCH_URL)

        try:
            data = json.loads(raw) if raw.strip() else {}
        except json.JSONDecodeError as exc:
            raise RuntimeError(f"RCSB search returned malformed JSON: {exc}") from exc

        hits: list[dict] = []
        for item in data.get("result_set", [])[: self.settings.rcsb_search_max_hits]:
            identifier = str(item.get("identifier", ""))
            pdb_id = identifier.split("_")[0].upper()
            if not is_pdb_id(pdb_id):
                continue
            entity_id = identifier.split("_")[1] if "_" in identifier else "1"
            hits.append(
                {
                    "pdb_id": pdb_id,
                    "entity_id": entity_id,
                    "score": float(item.get("score", 0.0)),
                    "accessions": self.entity_accessions(pdb_id, entity_id, refresh=refresh),
                }
            )
        return hits

    def entity_accessions(self, pdb_id: str, entity_id: str, *, refresh: bool = False) -> list[str]:
        """UniProt accessions cross-referenced by one polymer entity.

        A failure here is not fatal: without cross-references the mapping is
        simply reported at a lower confidence level rather than as an exact
        match, which is the conservative outcome.
        """
        key = self._cache_key(f"{pdb_id}_{entity_id}", kind="polymer")
        raw = self._cache_get(key, refresh)
        if raw is None:
            url = RCSB_POLYMER_URL.format(pdb_id=pdb_id, entity_id=entity_id)
            outcome = self.client.get(url)
            if not outcome.ok or outcome.response is None:
                return []
            raw = outcome.response.text
            self._cache_put(key, raw, f"{pdb_id}_{entity_id}", url)
        try:
            data = json.loads(raw)
        except json.JSONDecodeError:
            return []
        accessions: list[str] = []
        for container in data.get("rcsb_polymer_entity_container_identifiers", {}).get(
            "reference_sequence_identifiers", []
        ):
            if container.get("database_name") in ("UniProt", "UniProtKB"):
                accession = container.get("database_accession")
                if accession:
                    accessions.append(accession.upper())
        for aligned in data.get("rcsb_polymer_entity_align", []):
            if aligned.get("reference_database_name") in ("UniProt", "UniProtKB"):
                accession = aligned.get("reference_database_accession")
                if accession and accession.upper() not in accessions:
                    accessions.append(accession.upper())
        return accessions

    def entry_metadata(self, pdb_id: str, *, refresh: bool = False) -> dict:
        """Title, experimental method and resolution for an entry."""
        key = self._cache_key(pdb_id, kind="entry")
        raw = self._cache_get(key, refresh)
        if raw is None:
            url = RCSB_ENTRY_URL.format(pdb_id=pdb_id.upper())
            outcome = self.client.get(url)
            if not outcome.ok or outcome.response is None:
                return {}
            raw = outcome.response.text
            self._cache_put(key, raw, pdb_id, url)
        try:
            data = json.loads(raw)
        except json.JSONDecodeError:
            return {}
        methods = [m.get("method", "") for m in data.get("exptl", [])]
        resolution = (data.get("rcsb_entry_info") or {}).get("resolution_combined") or []
        return {
            "title": (data.get("struct") or {}).get("title", ""),
            "methods": methods,
            "resolution_a": resolution[0] if resolution else None,
        }

    def _search_payload(self, sequence: str) -> dict:
        return {
            "query": {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                    "evalue_cutoff": 1.0,
                    "identity_cutoff": self.settings.rcsb_search_identity_cutoff,
                    "sequence_type": "protein",
                    "value": sequence,
                },
            },
            "request_options": {
                "scoring_strategy": "sequence",
                "paginate": {"start": 0, "rows": self.settings.rcsb_search_max_hits},
            },
            "return_type": "polymer_entity",
        }

    # -- shared -------------------------------------------------------------

    def _finalise(
        self, result: StructureResult, pdb_text: str, pdb_id: str, url: str, from_cache: bool
    ) -> None:
        if not looks_like_pdb(pdb_text):
            result.status = RetrievalStatus.SERVICE_ERROR
            result.error = f"RCSB returned no coordinates for {pdb_id}."
            return
        parsed = parse_pdb(pdb_text)
        confidence, kind = mean_confidence(pdb_text, is_predicted=False)
        result.status = RetrievalStatus.OK
        result.identifier = pdb_id
        result.structure_url = url
        result.download_path = self._save(pdb_text, f"{pdb_id}.pdb")
        result.sequence = parsed.longest_chain.sequence if parsed.longest_chain else None
        result.confidence = confidence
        result.confidence_kind = kind
        result.experimental = True
        result.from_cache = from_cache


__all__ = ["RCSBProvider"]
