#!/usr/bin/env python3
# AMENDMENT-6
"""Experiment D: identity x coverage x accession, against INDEPENDENT truth.

Amendment 6 changes, found when the first corrected run was inspected:
  * Thresholds are read from bioseqinsight's own Settings (exact 99/95, high
    95/80, related 40/40). The previous version hard-coded related 30/50,
    which disagreed with the software and with nothing else.
  * The unrelated-pair noise floor is reported CONDITIONED ON COVERAGE.
    Identity is measured over the aligned region only, so a short subject
    chain gets a high chance identity at low coverage; the unconditioned
    percentile (p90 = 69%) therefore overstates how similar unrelated
    proteins look when a substantial part of the query is covered.

positives  = naive top hits that independent ground truth says ARE the entry
negatives  = naive top hits that it says are NOT, plus constructed negatives
             B1 (paralogue), B4 (unrelated), B5 (near-identical, other entry).
             B2 (right protein, partial structure) is reported separately;
             B3 (isoform) is excluded (its label is a construction assumption).
identity / coverage / accession_match are bioseqinsight's own measurements of
each hit (the system under test), not the labels.

Usage
-----
    python validation/scripts/04_threshold_analysis.py
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(VALIDATION_ROOT))
from lib import stats  # noqa: E402

POS = {"positive", "positive_single_source"}
NEG = {"negative", "negative_single_source"}
IDENTITY_BINS = list(range(0, 101, 10))
COVERAGE_BINS = list(range(0, 101, 10))
DEFAULTS = {"exact_identity": 99.0, "exact_coverage": 95.0, "high_identity": 95.0,
            "high_coverage": 80.0, "related_identity": 40.0, "related_coverage": 40.0}


def thresholds() -> dict:
    """bioseqinsight's own defaults when importable, else the documented ones."""
    try:
        repo_src = VALIDATION_ROOT.parent / "src"
        if repo_src.is_dir():
            sys.path.insert(0, str(repo_src))
        from bioseqinsight.config.settings import Settings
        s = Settings()
        return {k: float(getattr(s, f"{k}_threshold")) for k in DEFAULTS}
    except Exception:  # noqa: BLE001
        return dict(DEFAULTS)


def _f(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _tf(value) -> bool:
    return str(value).lower() == "true"


def load_points(baselines_dir: Path, adversarial_csv: Path) -> list[dict]:
    points: list[dict] = []
    path = baselines_dir / "rcsb_baseline_results.csv"
    if path.is_file():
        with open(path, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                identity, coverage = _f(row.get("method1_identity")), _f(row.get("method1_coverage"))
                label = row.get("method1_truth_label", "")
                if identity is None or coverage is None or (label not in POS and label not in NEG):
                    continue
                points.append({"identity": identity, "coverage": coverage, "positive": label in POS,
                               "source": "naive_hit", "accession_match": _tf(row.get("method1_accession_match")),
                               "accession": row.get("accession")})
    if adversarial_csv.is_file():
        with open(adversarial_csv, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                cls = row.get("class", "")
                if row.get("status") != "ok" or cls.startswith("B3"):
                    continue
                identity, coverage = _f(row.get("identity_percent")), _f(row.get("coverage_percent"))
                if identity is None or coverage is None:
                    continue
                points.append({"identity": identity, "coverage": coverage, "positive": False,
                               "source": cls, "accession_match": _tf(row.get("accession_match")),
                               "accession": row.get("query_accession")})
    return points


def _bin(value: float, edges: list[int]) -> int:
    for i in range(len(edges) - 1):
        if edges[i] <= value < edges[i + 1]:
            return i
    return len(edges) - 2


def _pct(sorted_values: list[float], q: float) -> float:
    return round(sorted_values[min(len(sorted_values) - 1, int(q * len(sorted_values)))], 2)


def analyse(points: list[dict], t: dict | None = None) -> dict:
    t = t or thresholds()
    exact = lambda p: p["identity"] >= t["exact_identity"] and p["coverage"] >= t["exact_coverage"]  # noqa: E731
    naive_pos = [p for p in points if p["source"] == "naive_hit" and p["positive"]]
    naive_neg = [p for p in points if p["source"] == "naive_hit" and not p["positive"]]
    adv_neg = [p for p in points if p["source"] in ("B1_wrong_paralogue", "B4_unrelated",
                                                    "B5_high_identity_wrong_entity")]
    fragments = [p for p in points if p["source"].startswith("B2")]

    def aucs(pos, neg):
        if not pos or not neg:
            return None
        return {"n_positive": len(pos), "n_negative": len(neg),
                "auc_identity": round(stats.auc([p["identity"] for p in pos], [p["identity"] for p in neg]), 4),
                "auc_coverage": round(stats.auc([p["coverage"] for p in pos], [p["coverage"] for p in neg]), 4),
                "auc_identity_x_coverage": round(stats.auc(
                    [p["identity"] * p["coverage"] / 100 for p in pos],
                    [p["identity"] * p["coverage"] / 100 for p in neg]), 4)}

    all_neg = naive_neg + adv_neg
    ic_pass = [p for p in all_neg if exact(p)]
    ablation = {
        "negatives_total": len(all_neg),
        "pass_identity_and_coverage_only": len(ic_pass),
        "pass_identity_coverage_and_accession": sum(1 for p in ic_pass if p["accession_match"]),
        "positives_total": len(naive_pos),
        "positives_pass_identity_and_coverage": sum(1 for p in naive_pos if exact(p)),
        "positives_pass_all_three": sum(1 for p in naive_pos if exact(p) and p["accession_match"]),
    }

    unrelated = [p for p in points if p["source"] == "B4_unrelated"]
    noise = None
    if unrelated:
        covered = [p for p in unrelated if p["coverage"] >= t["related_coverage"]]
        ids_all = sorted(p["identity"] for p in unrelated)
        ids_cov = sorted(p["identity"] for p in covered)
        reach = lambda i, c: sum(1 for p in unrelated if p["identity"] >= i and p["coverage"] >= c)  # noqa: E731
        noise = {
            "n_unrelated_pairs": len(unrelated),
            "identity_percentiles_all_pairs": {f"p{q}": _pct(ids_all, q / 100) for q in (50, 90, 95, 99)},
            "n_with_coverage_at_least_related_threshold": len(covered),
            "identity_percentiles_given_coverage_ge_related": (
                {f"p{q}": _pct(ids_cov, q / 100) for q in (50, 90, 95, 99)} if ids_cov else None),
            "high_identity_low_coverage_pairs": sum(
                1 for p in unrelated if p["identity"] >= 50 and p["coverage"] < t["related_coverage"]),
            "reach_M2_criteria": reach(t["related_identity"], t["related_coverage"]),
            "reach_M2_criteria_fraction": round(reach(t["related_identity"], t["related_coverage"]) / len(unrelated), 4),
            "reach_identity_50_and_related_coverage": reach(50, t["related_coverage"]),
            "reach_high_criteria": reach(t["high_identity"], t["high_coverage"]),
            "reach_exact_identity_and_coverage": reach(t["exact_identity"], t["exact_coverage"]),
            "note": ("Pairs here share no curated family. High identity at LOW coverage is a short-chain "
                     "artefact (identity is over the aligned region only); the conditioned percentiles are "
                     "the relevant noise floor for M2."),
        }

    grid_pos = [[0] * (len(COVERAGE_BINS) - 1) for _ in range(len(IDENTITY_BINS) - 1)]
    grid_neg = [[0] * (len(COVERAGE_BINS) - 1) for _ in range(len(IDENTITY_BINS) - 1)]
    for p in naive_pos + all_neg:
        (grid_pos if p["positive"] else grid_neg)[_bin(p["identity"], IDENTITY_BINS)][_bin(p["coverage"], COVERAGE_BINS)] += 1

    return {
        "thresholds_used": t,
        "naive_positive_vs_naive_negative": aucs(naive_pos, naive_neg),
        "naive_positive_vs_adversarial_negative": aucs(naive_pos, adv_neg),
        "ablation_at_M4_thresholds": ablation,
        "unrelated_pair_noise_floor": noise,
        "fragments_B2": {"n": len(fragments),
                         "median_identity": stats.median([p["identity"] for p in fragments]) if fragments else None,
                         "median_coverage": stats.median([p["coverage"] for p in fragments]) if fragments else None,
                         "accession_confirmed": sum(1 for p in fragments if p["accession_match"])},
        "interpretation": ("Coverage measures completeness, not whether the molecule is right: correct hits "
                           "include many fragments, so identity x coverage can score WORSE than identity alone. "
                           "The accession check, not coverage, removes the high-identity wrong-entry hits."),
        "grid_positive": grid_pos, "grid_negative": grid_neg,
        "identity_bins": IDENTITY_BINS, "coverage_bins": COVERAGE_BINS,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--baselines", default=str(VALIDATION_ROOT / "results" / "baselines"))
    parser.add_argument("--adversarial",
                        default=str(VALIDATION_ROOT / "results" / "adversarial" / "adversarial_results.csv"))
    parser.add_argument("--out", default=str(VALIDATION_ROOT / "results" / "threshold_analysis"))
    args = parser.parse_args(argv)
    points = load_points(Path(args.baselines), Path(args.adversarial))
    if not points:
        print("No points found. Run 02_run_baselines.py and 03_run_adversarial.py first.", file=sys.stderr)
        return 1
    result = analyse(points)
    out = Path(f"{args.out}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2))
    print(json.dumps({k: v for k, v in result.items() if not k.startswith("grid")}, indent=2))
    print(f"\nWrote {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
