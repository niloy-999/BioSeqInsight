# AMENDMENT-5
"""Analyse one candidate structure ("hit") against one query.

Used by Experiment C (naive baselines), D (threshold analysis) and E (case
studies) so that all three share one definition of what is measured:

  * ground truth   -- independent (lib/ground_truth.py: PDBe SIFTS + RCSB
                      polymer-entity references), NEVER bioseqinsight's own
                      alignment;
  * the measurement -- bioseqinsight's own validate_structure(), the system
                      under test, giving identity, coverage and an M-level.

Keeping those two apart is what prevents the circularity found in the first
full run (see PROTOCOL.md Amendment 5).
"""

from __future__ import annotations

from . import api, ground_truth

M_RANK = {"M0": 0, "M1": 1, "M2": 2, "M3": 3, "M4": 4}


def split_identifier(identifier: str) -> tuple[str, str]:
    """``"1UBQ_1"`` -> ``("1UBQ", "1")``; a bare PDB id gets entity ``"1"``."""
    pdb_id, _, entity_id = identifier.partition("_")
    return pdb_id.upper(), (entity_id or "1")


def hit_accessions(pdb_id: str, entity_id: str, use_cache: bool = True) -> list[str]:
    """UniProt accessions RCSB attaches to the hit entity (for human inspection)."""
    try:
        entity = api.rcsb_polymer_entity(pdb_id, entity_id, use_cache=use_cache)
    except api.FetchError:
        return []
    refs = (entity.get("rcsb_polymer_entity_container_identifiers", {})
            .get("reference_sequence_identifiers", []))
    return sorted({r["database_accession"] for r in refs
                   if r.get("database_name") == "UniProt" and r.get("database_accession")})


def analyse_hit(pdb_id: str, entity_id: str, accessions: list[str], sequence: str,
                settings, use_cache: bool = True, validate=None) -> dict:
    """Independent truth label + bioseqinsight's identity/coverage/M-level."""
    if validate is None:
        from bioseqinsight.structures.validation import validate_structure as validate

    truth = ground_truth.adjudicate_any(pdb_id, accessions, query_length=len(sequence),
                                        use_cache=use_cache)
    out: dict = {
        "pdb_id": pdb_id,
        "entity_id": entity_id,
        "truth_label": truth.label,
        "truth_sifts": truth.sifts_says_yes,
        "truth_rcsb": truth.rcsb_says_yes,
        "hit_accessions": ";".join(hit_accessions(pdb_id, entity_id, use_cache)),
    }
    try:
        pdb_text = api.fetch_rcsb_pdb_file(pdb_id, use_cache=use_cache)
    except api.FetchError as exc:
        out["validation_error"] = str(exc)
        return out

    best = None
    for accession in accessions:
        report = validate(pdb_text, query_sequence=sequence, query_accession=accession,
                          settings=settings)
        if best is None or M_RANK[report.mapping_level.value] > M_RANK[best.mapping_level.value]:
            best = report
    out["identity"] = best.identity_percent
    out["coverage"] = best.coverage_percent
    out["mapping"] = best.mapping_level.value
    out["accession_match"] = best.accession_match
    return out


__all__ = ["M_RANK", "split_identifier", "hit_accessions", "analyse_hit"]
