#!/usr/bin/env python3
# AMENDMENT-5
"""Experiment C (corrected): naive retrieval versus validated retrieval.

What was wrong before (PROTOCOL.md Amendment 5): Methods 1 and 2 were scored
as "accepted" exactly when their hit was ground-truth positive, and then
evaluated against that same ground truth -- so they scored a perfect 1.0 by
construction. A naive method accepts WHATEVER THE RESOURCE RETURNS.

Definitions used now (all per query protein, RCSB only; AlphaFold DB and ESM
Atlas are accession/sequence-keyed and are summarised separately):

  method1  naive      accept the top RCSB sequence-search hit (no identity cutoff)
  method2  naive      accept the top hit of a search carrying RCSB's own
                      identity cutoff (0.9)
  method1/2 validated the SAME hit, accepted only if bioseqinsight's validation
                      layer classifies it M4  -> isolates the value of the layer
  method3  end-to-end bioseqinsight's full retrieval (RCSB provider), accepted
                      only if the result is M4

Correctness of any accepted structure is judged by independent ground truth
(lib/ground_truth.py), never by bioseqinsight's alignment. ``has_structure``
(UniProt lists >=1 PDB entry) separates proteins for which a correct structure
exists from those for which none does -- for the latter, accepting anything
is a false association by definition.

Usage
-----
    python validation/scripts/02_run_baselines.py \\
        --dataset validation/data/dataset_accessions.csv --out validation/results/baselines
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(VALIDATION_ROOT))
from lib import api, ground_truth, hits  # noqa: E402
from lib.fasta import read_fasta  # noqa: E402

REPO_ROOT = VALIDATION_ROOT.parent
try:
    from bioseqinsight.config.settings import Settings
    from bioseqinsight.structures.manager import StructureManager
except ImportError:
    sys.path.insert(0, str(REPO_ROOT / "src"))
    from bioseqinsight.config.settings import Settings
    from bioseqinsight.structures.manager import StructureManager

RCSB_IDENTITY_CUTOFF = 0.9
CHECKPOINT_EVERY = 50


def search_top_hit(sequence: str, identity_cutoff: float | None, use_cache: bool):
    found = api.rcsb_sequence_search(sequence, identity_cutoff=identity_cutoff, use_cache=use_cache)
    if not found:
        return None
    return hits.split_identifier(found[0]["identifier"])


def run_naive_methods(accession: str, sequence: str, settings, use_cache: bool, validate=None) -> dict:
    row: dict = {}
    for label, cutoff in (("method1", None), ("method2", RCSB_IDENTITY_CUTOFF)):
        try:
            hit = search_top_hit(sequence, cutoff, use_cache)
        except api.FetchError as exc:
            row[f"{label}_hit"] = None
            row[f"{label}_error"] = str(exc)
            continue
        if hit is None:
            row[f"{label}_hit"] = False
            continue
        pdb_id, entity_id = hit
        analysis = hits.analyse_hit(pdb_id, entity_id, [accession], sequence, settings,
                                    use_cache=use_cache, validate=validate)
        row[f"{label}_hit"] = True            # a naive method accepts what it is given
        for key, value in analysis.items():
            row[f"{label}_{key}"] = value
        row[f"{label}_validated_accept"] = analysis.get("mapping") == "M4"
    return row


def run_method3(manager, accession: str, sequence: str, resource: str, use_cache: bool) -> dict:
    outcome = manager.retrieve(accession, providers=[resource], stop_on_exact=False)
    best = outcome.best
    if best is None:
        return {"method3_status": "no_result", "method3_mapping": "M0", "method3_accepted": False}
    result = {
        "method3_status": "ok",
        "method3_mapping": best.mapping_level.value,
        "method3_identifier": best.identifier,
        "method3_accepted": best.mapping_level.value == "M4",
    }
    if best.validation:
        result["method3_identity"] = best.validation.identity_percent
        result["method3_coverage"] = best.validation.coverage_percent
        result["method3_accession_match"] = best.validation.accession_match
    if resource == "rcsb" and best.identifier:
        pdb_id = best.identifier.split("_")[0]
        truth = ground_truth.adjudicate_any(pdb_id, [accession], query_length=len(sequence),
                                            use_cache=use_cache)
        result["method3_truth_label"] = truth.label
    return result


def write_csv(path: Path, rows: list[dict]) -> None:
    fieldnames = ["accession"] + sorted({k for r in rows for k in r} - {"accession"})
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--fasta")
    parser.add_argument("--out", default=str(VALIDATION_ROOT / "results" / "baselines"))
    parser.add_argument("--limit", type=int)
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--config")
    args = parser.parse_args(argv)
    use_cache = not args.no_cache

    fasta_path = Path(args.fasta) if args.fasta else Path(args.dataset).with_name(
        Path(args.dataset).name.replace("_accessions.csv", "_sequences.fasta"))
    sequences = {acc: seq for acc, _d, seq in read_fasta(fasta_path)}
    with open(args.dataset, newline="", encoding="utf-8") as fh:
        candidates = list(csv.DictReader(fh))
    if args.limit:
        candidates = candidates[: args.limit]

    settings = Settings.load(path=args.config)
    manager = StructureManager(settings)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    rcsb_rows: list[dict] = []
    af_rows: list[dict] = []
    started = time.time()

    def flush() -> None:
        if rcsb_rows:
            write_csv(out_dir / "rcsb_baseline_results.csv", rcsb_rows)
        if af_rows:
            write_csv(out_dir / "alphafold_baseline_results.csv", af_rows)

    for i, cand in enumerate(candidates, 1):
        accession = cand["accession"]
        sequence = sequences.get(accession)
        print(f"[{i}/{len(candidates)}] {accession}", file=sys.stderr)
        base = {"accession": accession, "length": cand.get("length"),
                "has_structure": bool(cand.get("pdb_ids"))}

        if sequence:
            row = dict(base)
            try:
                row.update(run_naive_methods(accession, sequence, settings, use_cache))
                row.update(run_method3(manager, accession, sequence, "rcsb", use_cache))
            except Exception as exc:  # noqa: BLE001 - one candidate must not abort a long run
                row["error"] = f"{type(exc).__name__}: {exc}"
                print(f"  ! {accession}: {row['error']}", file=sys.stderr)
            rcsb_rows.append(row)

            af_row = dict(base)
            try:
                af_row.update(run_method3(manager, accession, sequence, "alphafold", use_cache))
            except Exception as exc:  # noqa: BLE001
                af_row["error"] = f"{type(exc).__name__}: {exc}"
                print(f"  ! {accession} (alphafold): {af_row['error']}", file=sys.stderr)
            af_rows.append(af_row)

        if i % CHECKPOINT_EVERY == 0:
            flush()

    flush()
    meta = {"n_candidates": len(candidates), "elapsed_s": round(time.time() - started, 1),
            "rcsb_identity_cutoff_method2": RCSB_IDENTITY_CUTOFF,
            "definition": "naive methods accept whatever the resource returns (Amendment 5)"}
    (out_dir / "run_metadata.json").write_text(json.dumps(meta, indent=2))
    print(json.dumps(meta, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
