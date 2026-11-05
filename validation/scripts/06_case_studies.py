#!/usr/bin/env python3
# AMENDMENT-5
"""Experiment E (corrected): adjudicate the v1.0 25-protein benchmark.

The 25 queries and the "mapped ..." notes are quoted verbatim from Table 3 of
the v1.0 manuscript (archived at tag v1.0-archive). v1.0's "AF pLDDT" column
is AlphaFold's confidence in its own model for whichever accession the
retrieval landed on; it is NOT identity to the query.

What was wrong before (PROTOCOL.md Amendment 5): the first version asked
bioseqinsight to retrieve each ACCESSION. That is largely answered by
AlphaFold DB, which is accession-keyed and correct by construction, so all 7
discordant queries trivially returned M4. That does not say whether v1.0's
discordant HITS were right.

Now, for each query this script reproduces what a sequence-based retrieval
returns TODAY (database contents change, so this is not v1.0's historical
2026-08-28 result) and judges it:

  * the top RCSB sequence-search hit for the query sequence (no cutoff);
  * its UniProt accession(s) as curated by RCSB (e.g. an IDE entry for insulin);
  * independent ground truth: does that ENTRY contain the queried accession
    (any chain), per PDBe SIFTS and RCSB entity references;
  * bioseqinsight's identity / coverage / M-level for that same hit;
  * for reference, bioseqinsight's accession-keyed retrieval from RCSB alone
    and from AlphaFold DB alone.

The 76-residue ubiquitin sequence is encoded identically by four genes
(P0CG48, P0CG47, P62979, P62987); an entry mapped to any of them counts.

Usage
-----
    python validation/scripts/06_case_studies.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(VALIDATION_ROOT))
from lib import api, hits  # noqa: E402

REPO_ROOT = VALIDATION_ROOT.parent
try:
    from bioseqinsight.config.settings import Settings
    from bioseqinsight.structures.manager import StructureManager
except ImportError:
    sys.path.insert(0, str(REPO_ROOT / "src"))
    from bioseqinsight.config.settings import Settings
    from bioseqinsight.structures.manager import StructureManager

UBIQUITIN_SEQUENCE = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
UBIQUITIN_ACCESSIONS = ["P0CG48", "P0CG47", "P62979", "P62987"]

V1_TABLE_3 = [
    {"query_name": "BPTI", "accession": "P00974", "length": 100, "esm_status": "504", "af_plddt": 97.80, "hit_note": "mapped Q52V24"},
    {"query_name": "SUMO1", "accession": "P63165", "length": 101, "esm_status": "504", "af_plddt": 78.28, "hit_note": None},
    {"query_name": "GroES", "accession": "P0A6F9", "length": 97, "esm_status": "504", "af_plddt": 90.80, "hit_note": "mapped GroEL"},
    {"query_name": "lysozyme (H, human)", "accession": "P61626", "length": 148, "esm_status": "ok", "af_plddt": 94.07, "hit_note": None},
    {"query_name": "lysozyme (hen)", "accession": "P00698", "length": 147, "esm_status": "ok", "af_plddt": 93.89, "hit_note": None},
    {"query_name": "T4 lysozyme", "accession": "P00720", "length": 164, "esm_status": "ok", "af_plddt": None, "hit_note": "AFDB 404"},
    {"query_name": "beta-2 microglobulin", "accession": "P61769", "length": 119, "esm_status": "ok", "af_plddt": 86.82, "hit_note": "mapped Q30597"},
    {"query_name": "insulin (preproinsulin)", "accession": "P01308", "length": 110, "esm_status": "ok", "af_plddt": 93.99, "hit_note": "mapped insulin-degrading enzyme (IDE)"},
    {"query_name": "RNase A", "accession": "P61823", "length": 150, "esm_status": "ok", "af_plddt": 94.04, "hit_note": None},
    {"query_name": "myoglobin", "accession": "P02144", "length": 154, "esm_status": "ok", "af_plddt": 97.18, "hit_note": None},
    {"query_name": "haemoglobin alpha", "accession": "P69905", "length": 142, "esm_status": "ok", "af_plddt": 98.06, "hit_note": None},
    {"query_name": "haemoglobin beta", "accession": "P68871", "length": 147, "esm_status": "ok", "af_plddt": 98.06, "hit_note": "mapped Hb-alpha"},
    {"query_name": "transthyretin", "accession": "P02766", "length": 147, "esm_status": "504", "af_plddt": 87.99, "hit_note": None},
    {"query_name": "CheY", "accession": "P0AE67", "length": 129, "esm_status": "504", "af_plddt": 96.99, "hit_note": None},
    {"query_name": "cyclophilin A", "accession": "P62937", "length": 165, "esm_status": "504", "af_plddt": 98.06, "hit_note": None},
    {"query_name": "H-Ras", "accession": "P01112", "length": 189, "esm_status": "504", "af_plddt": 91.95, "hit_note": None},
    {"query_name": "GFP", "accession": "P42212", "length": 238, "esm_status": "ok", "af_plddt": 96.64, "hit_note": None},
    {"query_name": "carbonic anhydrase II", "accession": "P00918", "length": 260, "esm_status": "ok", "af_plddt": 97.37, "hit_note": None},
    {"query_name": "IL-1 beta", "accession": "P01584", "length": 269, "esm_status": "ok", "af_plddt": 76.23, "hit_note": None},
    {"query_name": "CDK2", "accession": "P24941", "length": 298, "esm_status": "ok", "af_plddt": 88.46, "hit_note": "+ PDB 1AQ1"},
    {"query_name": "TNF-alpha", "accession": "P01375", "length": 233, "esm_status": "ok", "af_plddt": 84.56, "hit_note": None},
    {"query_name": "FGF2", "accession": "P09038", "length": 288, "esm_status": "504", "af_plddt": 70.27, "hit_note": None},
    {"query_name": "p53", "accession": "P04637", "length": 393, "esm_status": "ok", "af_plddt": 75.99, "hit_note": "mapped RPB1"},
    {"query_name": "MBP", "accession": "P0AEX9", "length": 396, "esm_status": "504", "af_plddt": 94.28, "hit_note": None},
    {"query_name": "ubiquitin (76 aa)", "accession": None, "length": 76, "esm_status": "ok", "af_plddt": 88.65, "hit_note": "mapped P0CG48"},
]


def is_discordant(entry: dict) -> bool:
    """v1.0's own "mapped <accession>" flag; other notes (e.g. "+ PDB 1AQ1",
    "AFDB 404") are not an accession discordance."""
    return (entry["hit_note"] or "").strip().lower().startswith("mapped ")


def query_for(entry: dict) -> tuple[list[str], str]:
    if entry["accession"] is None:
        return UBIQUITIN_ACCESSIONS, UBIQUITIN_SEQUENCE
    sequence = api.uniprot_entry(entry["accession"])["sequence"]["value"]
    return [entry["accession"]], sequence


def adjudicate_one(manager, settings, entry: dict, use_cache: bool = True, validate=None) -> dict:
    row = {**entry, "v1_flagged_discordant": is_discordant(entry)}
    accessions, sequence = query_for(entry)
    row["query_length_2026"] = len(sequence)

    found = api.rcsb_sequence_search(sequence, identity_cutoff=None, use_cache=use_cache)
    if not found:
        row["top_hit"] = "none"
    else:
        pdb_id, entity_id = hits.split_identifier(found[0]["identifier"])
        analysis = hits.analyse_hit(pdb_id, entity_id, accessions, sequence, settings,
                                    use_cache=use_cache, validate=validate)
        row["top_hit"] = "found"
        for key, value in analysis.items():
            row[f"top_hit_{key}"] = value

    primary = accessions[0]
    for resource in ("rcsb", "alphafold"):
        outcome = manager.retrieve(primary, providers=[resource], stop_on_exact=False)
        best = outcome.best
        row[f"accession_retrieval_{resource}_mapping"] = best.mapping_level.value if best else "M0"
        row[f"accession_retrieval_{resource}_id"] = best.identifier if best else ""
    return row


def summarise(rows: list[dict]) -> dict:
    discordant = [r for r in rows if r.get("v1_flagged_discordant")]

    def count(group, key, value):
        return sum(1 for r in group if r.get(key) == value)

    def block(group):
        truth = lambda r: r.get("top_hit_truth_label", "")  # noqa: E731
        return {
            "n": len(group),
            "top_hit_found": count(group, "top_hit", "found"),
            "top_hit_entry_contains_query_accession": sum(1 for r in group if truth(r) in ("positive", "positive_single_source")),
            "top_hit_entry_does_not_contain_query_accession": sum(1 for r in group if truth(r) in ("negative", "negative_single_source")),
            "top_hit_ambiguous_or_unknown": sum(1 for r in group if truth(r) in ("ambiguous", "unknown")),
            "top_hit_validated_M4": count(group, "top_hit_mapping", "M4"),
            "top_hit_validated_M3": count(group, "top_hit_mapping", "M3"),
            "top_hit_validated_below_M3": sum(1 for r in group if r.get("top_hit_mapping") in ("M0", "M1", "M2")),
        }

    return {"v1_discordant_queries": block(discordant),
            "v1_other_queries": block([r for r in rows if not r.get("v1_flagged_discordant")]),
            "interpretation_note": (
                "'Contains the query accession' is an accession-level claim. A top hit that is an "
                "identical or near-identical protein filed under a different accession, or a complex "
                "whose queried chain is not catalogued under this accession, counts as NOT containing "
                "it; inspect each such row (top_hit_hit_accessions) before describing it as an error.")}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(VALIDATION_ROOT / "results" / "case_studies"))
    parser.add_argument("--config")
    parser.add_argument("--no-cache", action="store_true")
    args = parser.parse_args(argv)

    settings = Settings.load(path=args.config)
    manager = StructureManager(settings)
    rows = []
    for i, entry in enumerate(V1_TABLE_3, 1):
        print(f"[{i}/{len(V1_TABLE_3)}] {entry['query_name']} ({entry['accession']})", file=sys.stderr)
        try:
            rows.append(adjudicate_one(manager, settings, entry, use_cache=not args.no_cache))
        except Exception as exc:  # noqa: BLE001
            print(f"  ! error: {exc}", file=sys.stderr)
            rows.append({**entry, "v1_flagged_discordant": is_discordant(entry), "error": f"{type(exc).__name__}: {exc}"})

    out_csv = Path(f"{args.out}.csv")
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({k for r in rows for k in r})
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    summary = summarise(rows)
    Path(f"{args.out}_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {out_csv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
