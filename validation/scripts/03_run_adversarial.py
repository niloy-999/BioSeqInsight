#!/usr/bin/env python3
"""Experiment B execution: run every adversarial case through validation
directly (PROTOCOL.md Section 3), bypassing bioseqinsight's own search.

Each adversarial case names a query and a *forced* structure the search
step would never have found on its own (by construction). This script
fetches that forced structure's real coordinate file independently (never
through bioseqinsight's own retrieval/search path, since these pairs are
deliberately not what a search would return) and calls bioseqinsight's own
``validate_structure`` function directly on the (query, forced structure)
pair -- which is exactly the function under test, used the way the real
pipeline uses it once a candidate structure is in hand.

By construction, NONE of these cases should classify as M4; the outcome
recorded for each is whether the system correctly stayed at or below its
``expected_max_m_level`` (never over-confident) rather than whether it hit
an exact label (adversarial construction estimates, e.g. the k-mer identity
in B1/B5, are approximate by design -- see 01_build_adversarial.py -- so the
precise M-level is informative but the "never M4" check is the real test).

Usage
-----
    python validation/scripts/03_run_adversarial.py \\
        --cases validation/data/adversarial_cases.csv \\
        --dataset validation/data/dataset_accessions.csv \\
        --out validation/results/adversarial
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(VALIDATION_ROOT))
from lib import api  # noqa: E402
from lib.fasta import read_fasta  # noqa: E402

REPO_ROOT = VALIDATION_ROOT.parent
try:
    from bioseqinsight.config.settings import Settings
    from bioseqinsight.structures.validation import validate_structure
except ImportError:
    sys.path.insert(0, str(REPO_ROOT / "src"))
    from bioseqinsight.config.settings import Settings
    from bioseqinsight.structures.validation import validate_structure

M_LEVEL_RANK = {"M0": 0, "M1": 1, "M2": 2, "M3": 3, "M4": 4}


def resolve_forced_pdb_id(case: dict, dataset_by_accession: dict[str, dict],
                           sequences: dict[str, str], use_cache: bool) -> str | None:
    """Find a concrete PDB ID for whatever the case's forced side specifies."""
    if case.get("forced_pdb_id"):
        return case["forced_pdb_id"]

    if case.get("class") == "B2_fragment":  # Amendment 4: these were silently skipped before
        row = dataset_by_accession.get(case["query_accession"])
        if row and row.get("pdb_ids"):
            return row["pdb_ids"].split(";")[0]
        return None

    if case.get("forced_structure_isoform_id"):
        iso_id = case["forced_structure_isoform_id"]
        try:
            iso_seq = api.uniprot_isoform_sequence(iso_id, use_cache=use_cache)
            hits = api.rcsb_sequence_search(iso_seq, identity_cutoff=0.95, use_cache=use_cache)
        except api.FetchError:
            return None
        if not hits:
            return None
        return hits[0]["identifier"].split("_")[0]

    forced_acc = case.get("forced_structure_accession")
    if forced_acc:
        row = dataset_by_accession.get(forced_acc)
        if row and row.get("pdb_ids"):
            return row["pdb_ids"].split(";")[0]
        seq = sequences.get(forced_acc)
        if seq:
            try:
                hits = api.rcsb_sequence_search(seq, identity_cutoff=None, use_cache=use_cache)
            except api.FetchError:
                return None
            if hits:
                return hits[0]["identifier"].split("_")[0]
    return None


def run_one(case: dict, settings, dataset_by_accession: dict, sequences: dict,
            use_cache: bool) -> dict:
    query_acc = case["query_accession"]
    query_seq = sequences.get(query_acc)
    row = {**case}
    if not query_seq:
        row["status"] = "no_query_sequence"
        return row

    pdb_id = resolve_forced_pdb_id(case, dataset_by_accession, sequences, use_cache)
    if not pdb_id:
        row["status"] = "no_forced_structure_resolved"
        return row
    row["resolved_pdb_id"] = pdb_id

    try:
        pdb_text = api.fetch_rcsb_pdb_file(pdb_id, use_cache=use_cache)
    except api.FetchError as exc:
        row["status"] = f"fetch_error: {exc}"
        return row

    report = validate_structure(
        pdb_text, query_sequence=query_seq, query_accession=query_acc, settings=settings,
    )
    row["status"] = "ok"
    row["mapping_level"] = report.mapping_level.value
    row["identity_percent"] = report.identity_percent
    row["coverage_percent"] = report.coverage_percent
    row["accession_match"] = report.accession_match
    row["subject_length"] = getattr(report, "subject_length", None)  # Amendment 4
    row["partial_structure"] = getattr(report, "partial_structure", None)  # Amendment 7
    row["query_length"] = getattr(report, "query_length", None)

    expected_max = case.get("expected_max_m_level", "M3")
    row["stayed_at_or_below_expected"] = (
        M_LEVEL_RANK.get(report.mapping_level.value, 0) <= M_LEVEL_RANK.get(expected_max, 3)
    )
    row["incorrectly_reached_m4"] = report.mapping_level.value == "M4"
    return row


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cases", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--fasta")
    parser.add_argument("--out", default=str(VALIDATION_ROOT / "results" / "adversarial"))
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--config")
    parser.add_argument("--limit", type=int)
    args = parser.parse_args(argv)
    use_cache = not args.no_cache

    with open(args.cases, newline="", encoding="utf-8") as fh:
        cases = list(csv.DictReader(fh))
    if args.limit:
        cases = cases[: args.limit]

    with open(args.dataset, newline="", encoding="utf-8") as fh:
        dataset_rows = list(csv.DictReader(fh))
    dataset_by_accession = {r["accession"]: r for r in dataset_rows}

    fasta_path = Path(args.fasta) if args.fasta else Path(args.dataset).with_name(
        Path(args.dataset).name.replace("_accessions.csv", "_sequences.fasta")
    )
    sequences = {acc: seq for acc, _desc, seq in read_fasta(fasta_path)}
    # B3 queries the canonical isoform but adversarial cases may also need
    # forced-isoform sequences resolved for identity reporting, not just PDB lookup.

    settings = Settings.load(path=args.config)

    rows = []
    for i, case in enumerate(cases, 1):
        print(f"[{i}/{len(cases)}] {case['class']}: {case['query_accession']}", file=sys.stderr)
        try:
            rows.append(run_one(case, settings, dataset_by_accession, sequences, use_cache))
        except Exception as exc:  # noqa: BLE001
            print(f"  ! error: {exc}", file=sys.stderr)
            rows.append({**case, "status": f"error: {exc}"})

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "adversarial_results.csv"
    fieldnames = sorted({k for row in rows for k in row})
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    ok_rows = [r for r in rows if r.get("status") == "ok"]
    by_class: dict[str, dict] = {}
    for row in ok_rows:
        cls = row["class"]
        bucket = by_class.setdefault(cls, {"n": 0, "stayed_within_expected": 0, "false_m4": 0})
        bucket["n"] += 1
        bucket["stayed_within_expected"] += int(row["stayed_at_or_below_expected"])
        bucket["false_m4"] += int(row["incorrectly_reached_m4"])

    summary = {
        "total_cases": len(rows), "resolved_and_scored": len(ok_rows),
        "unresolved": len(rows) - len(ok_rows),
        "by_class": by_class,
        "overall_false_m4_rate": (
            round(sum(b["false_m4"] for b in by_class.values()) / len(ok_rows), 4) if ok_rows else None
        ),
    }
    (out_dir / "adversarial_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {out_csv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
