# AMENDMENT-5
"""Independent ground truth: does PDB entry X correspond to UniProt accession Y?

PROTOCOL.md Section 5: ground truth must not depend on bioseqinsight's own
alignment engine. Two independent curated sources are cross-checked:

  1. PDBe SIFTS UniProt mapping (computed by EBI)
  2. RCSB's own polymer-entity reference-sequence identifiers

A case is POSITIVE only if both agree. Disagreement is "ambiguous", never
silently resolved either way. "Positive" means an ACCESSION-level claim: the
entry contains a chain that the curated databases map to that exact UniProt
accession. A structure of an identical or near-identical protein filed under
a different accession is NEGATIVE here by design (the question asked is
"is this the structure of this entry", which is what M4 asserts).

Amendment 5: adjudicate_entry()/adjudicate_any() check ALL polymer entities
of an entry. The earlier entity-specific check examined entity "1" only, so a
correct entry in which the queried protein is entity 2 (e.g. haemoglobin
beta) could be mislabelled.
"""

from __future__ import annotations

from dataclasses import dataclass

from . import api

POSITIVE = ("positive", "positive_single_source")
_RANK = {"positive": 5, "positive_single_source": 4, "ambiguous": 3,
         "unknown": 2, "negative_single_source": 1, "negative": 0}


@dataclass
class GroundTruthLabel:
    pdb_id: str
    accession: str
    sifts_says_yes: bool | None
    rcsb_says_yes: bool | None
    coverage_fraction: float | None
    coverage_is_fraction: bool
    label: str

    def to_row(self) -> dict:
        return {
            "pdb_id": self.pdb_id, "accession": self.accession,
            "sifts_says_yes": self.sifts_says_yes, "rcsb_says_yes": self.rcsb_says_yes,
            "coverage_fraction": self.coverage_fraction,
            "coverage_is_fraction": self.coverage_is_fraction, "label": self.label,
        }


def _sifts_agrees(pdb_id: str, accession: str, query_length: int | None = None,
                  use_cache: bool = True) -> tuple[bool | None, float | None, bool]:
    """PDB-ID-keyed PDBe lookup: ``/mappings/uniprot/{pdb_id}``.

    Returns ``(says_yes, coverage_or_span, coverage_is_fraction)``; ``None``
    means "unknown" (lookup failed or the response had an unexpected shape).
    """
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    try:
        body, _meta = api.cached_get_json(url, namespace="sifts_pdb_mapping",
                                          limiter=api._sifts_limiter, use_cache=use_cache)
    except api.FetchError as exc:
        if exc.status == 404:
            return False, None, False
        return None, None, False
    if not isinstance(body, dict):
        return None, None, False
    block = body.get(pdb_id.lower(), {})
    if not isinstance(block, dict):
        return None, None, False
    block = block.get("UniProt", {})
    if not isinstance(block, dict) or accession not in block:
        return False, None, False
    best_span = 0
    for mapping in block[accession].get("mappings", []):
        start = mapping.get("unp_start") or (mapping.get("start") or {}).get("author_residue_number")
        end = mapping.get("unp_end") or (mapping.get("end") or {}).get("author_residue_number")
        if start is not None and end is not None:
            best_span = max(best_span, end - start + 1)
    if not best_span:
        return True, None, False
    if query_length:
        return True, best_span / query_length, True
    return True, float(best_span), False


def _rcsb_agrees(pdb_id: str, entity_id: str, accession: str, use_cache: bool = True) -> bool | None:
    try:
        entity = api.rcsb_polymer_entity(pdb_id, entity_id, use_cache=use_cache)
    except api.FetchError:
        return None
    refs = (entity.get("rcsb_polymer_entity_container_identifiers", {})
            .get("reference_sequence_identifiers", []))
    accessions = {r.get("database_accession") for r in refs if r.get("database_name") == "UniProt"}
    if not accessions:
        return None
    return accession in accessions


def _rcsb_entry_agrees(pdb_id: str, accession: str, use_cache: bool = True) -> bool | None:
    """True if ANY polymer entity references the accession; False if entities
    carry UniProt references but none is this accession; None if undecidable."""
    try:
        entry = api.rcsb_entry(pdb_id, use_cache=use_cache)
    except api.FetchError:
        return None
    entity_ids = (entry.get("rcsb_entry_container_identifiers", {}) or {}).get("polymer_entity_ids") or ["1"]
    saw_refs = False
    for entity_id in entity_ids:
        verdict = _rcsb_agrees(pdb_id, str(entity_id), accession, use_cache=use_cache)
        if verdict is True:
            return True
        if verdict is False:
            saw_refs = True
    return False if saw_refs else None


def _combine(sifts_yes: bool | None, rcsb_yes: bool | None) -> str:
    if sifts_yes is None and rcsb_yes is None:
        return "unknown"
    if sifts_yes is True and rcsb_yes is True:
        return "positive"
    if sifts_yes is False and rcsb_yes is False:
        return "negative"
    if sifts_yes is None or rcsb_yes is None:
        return "positive_single_source" if (sifts_yes or rcsb_yes) else "negative_single_source"
    return "ambiguous"


def adjudicate(pdb_id: str, entity_id: str, accession: str, query_length: int | None = None,
               use_cache: bool = True) -> GroundTruthLabel:
    """Entity-specific adjudication (kept for compatibility)."""
    sifts_yes, coverage, is_fraction = _sifts_agrees(pdb_id, accession, query_length, use_cache)
    rcsb_yes = _rcsb_agrees(pdb_id, entity_id, accession, use_cache)
    return GroundTruthLabel(pdb_id, accession, sifts_yes, rcsb_yes, coverage, is_fraction,
                            _combine(sifts_yes, rcsb_yes))


def adjudicate_entry(pdb_id: str, accession: str, query_length: int | None = None,
                     use_cache: bool = True) -> GroundTruthLabel:
    """Does ANY chain/entity of this entry correspond to the accession?"""
    sifts_yes, coverage, is_fraction = _sifts_agrees(pdb_id, accession, query_length, use_cache)
    rcsb_yes = _rcsb_entry_agrees(pdb_id, accession, use_cache)
    return GroundTruthLabel(pdb_id, accession, sifts_yes, rcsb_yes, coverage, is_fraction,
                            _combine(sifts_yes, rcsb_yes))


def adjudicate_any(pdb_id: str, accessions: list[str], query_length: int | None = None,
                   use_cache: bool = True) -> GroundTruthLabel:
    """Best label over several acceptable accessions (e.g. identical ubiquitin
    sequences encoded by four genes)."""
    best: GroundTruthLabel | None = None
    for accession in accessions:
        label = adjudicate_entry(pdb_id, accession, query_length, use_cache)
        if best is None or _RANK[label.label] > _RANK[best.label]:
            best = label
    assert best is not None, "adjudicate_any needs at least one accession"
    return best


__all__ = ["GroundTruthLabel", "POSITIVE", "adjudicate", "adjudicate_entry", "adjudicate_any"]
