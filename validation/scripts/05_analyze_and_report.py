#!/usr/bin/env python3
# AMENDMENT-6
"""Analyse Experiment C (baseline comparison), two operating points.

Unit: one query protein. Each method either accepts a structure or returns
nothing. An accepted structure is judged by INDEPENDENT ground truth (PDBe
SIFTS + RCSB entity references), at ACCESSION level: correct (TP), wrong
(FP = a false association), or unadjudicable. Strata are never blended:
has_structure (a correct answer exists) / no_structure (accepting anything is
a false association by definition) / all.

Methods
  method1_naive / method2_naive      accept whatever the resource returns
  *_validated_strict                 the SAME hit, accepted only if M4
  *_validated_permissive             the SAME hit, accepted if M3 or M4
  method3_strict / method3_permissive bioseqinsight's own retrieval, M4 / M3+

Operating points (Amendment 6): STRICT (M4) asserts an exact match; PERMISSIVE
(M3+) also accepts a high-identity, high-coverage hit whose accession is not
confirmed, which can be the same sequence filed under another accession.

Because "false" means "not this accession", every false association is also
placed in a TIER by the hit's own measured identity/coverage, so a reader can
see how many are genuinely different proteins:
  same_sequence   identity >= 99 and coverage >= 95
  near_identical  identity >= 90 (otherwise)
  different       identity <  90
  unmeasured      no coordinate file could be validated

Also reported: why correct naive hits were NOT accepted by the strict gate
(the cost of the gate), and why hits were unmeasured.

Usage
-----
    python validation/scripts/05_analyze_and_report.py
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
M3_PLUS = ("M3", "M4")

METHODS = {
    "method1_naive": "top RCSB sequence-search hit, accepted unconditionally",
    "method1_validated_strict": "the same hit, accepted only if M4",
    "method1_validated_permissive": "the same hit, accepted if M3 or M4",
    "method2_naive": "top hit with RCSB's own 0.9 identity cutoff, accepted unconditionally",
    "method2_validated_strict": "the same hit, accepted only if M4",
    "method2_validated_permissive": "the same hit, accepted if M3 or M4",
    "method3_strict": "bioseqinsight retrieval, accepted only if M4",
    "method3_permissive": "bioseqinsight retrieval, accepted if M3 or M4",
}
TIERS = ("same_sequence", "near_identical", "different", "unmeasured")


def _tf(value) -> bool:
    return str(value).lower() == "true"


def _f(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def prefix(method: str) -> str:
    return method[:7]  # method1 / method2 / method3


def decision(row: dict, method: str) -> tuple[bool, str]:
    """(accepted?, ground-truth label of the accepted structure)."""
    p = prefix(method)
    if p == "method3":
        mapping = row.get("method3_mapping", "")
        accepted = mapping == "M4" if method.endswith("strict") else mapping in M3_PLUS
        return accepted, row.get("method3_truth_label", "")
    hit = _tf(row.get(f"{p}_hit"))
    label = row.get(f"{p}_truth_label", "")
    if method.endswith("naive"):
        return hit, label
    mapping = row.get(f"{p}_mapping", "")
    accepted = mapping == "M4" if method.endswith("strict") else mapping in M3_PLUS
    return hit and accepted, label


def outcome(row: dict, method: str) -> str:
    accepted, label = decision(row, method)
    if not accepted:
        return "none"
    if label in POS:
        return "tp"
    if label in NEG:
        return "fp"
    return "unadjudicable"


def tier(row: dict, method: str) -> str:
    p = prefix(method)
    identity, coverage = _f(row.get(f"{p}_identity")), _f(row.get(f"{p}_coverage"))
    if identity is None or coverage is None:
        return "unmeasured"
    if identity >= 99 and coverage >= 95:
        return "same_sequence"
    if identity >= 90:
        return "near_identical"
    return "different"


def _stratum(rows: list[dict], stratum: str) -> list[dict]:
    if stratum == "all":
        return rows
    return [r for r in rows if _tf(r.get("has_structure")) == (stratum == "has_structure")]


def summarise(rows: list[dict], method: str, stratum: str) -> dict:
    sub = _stratum(rows, stratum)
    counts = {"tp": 0, "fp": 0, "unadjudicable": 0, "none": 0}
    tiers = {t: 0 for t in TIERS}
    for r in sub:
        o = outcome(r, method)
        counts[o] += 1
        if o == "fp":
            tiers[tier(r, method)] += 1
    n, tp, fp = len(sub), counts["tp"], counts["fp"]
    far, far_lo, far_hi = stats.wilson(fp, tp + fp)
    result = {
        "n_queries": n, "accepted": tp + fp + counts["unadjudicable"],
        "true_associations": tp, "false_associations": fp,
        "false_association_tiers": tiers,
        "unadjudicable_accepted": counts["unadjudicable"], "returned_nothing": counts["none"],
        "false_association_rate": None if tp + fp == 0 else round(far, 4),
        "false_association_rate_95ci": None if tp + fp == 0 else [round(far_lo, 4), round(far_hi, 4)],
        "false_associations_per_1000_queries": round(1000 * fp / n, 1) if n else None,
    }
    if stratum == "has_structure" and n:
        s, s_lo, s_hi = stats.wilson(tp, n)
        result["sensitivity_correct_structure_returned"] = round(s, 4)
        result["sensitivity_95ci"] = [round(s_lo, 4), round(s_hi, 4)]
    return result


def mcnemar(rows: list[dict], a: str, b: str) -> dict:
    only_a = sum(1 for r in rows if outcome(r, a) == "fp" and outcome(r, b) != "fp")
    only_b = sum(1 for r in rows if outcome(r, b) == "fp" and outcome(r, a) != "fp")
    return {"queries": len(rows), "false_assoc_only_in_first": only_a,
            "false_assoc_only_in_second": only_b, "p_value": stats.mcnemar_exact(only_a, only_b)}


def loss_reason(row: dict, p: str = "method1") -> str:
    """Why a CORRECT naive hit was or was not accepted by the strict gate."""
    mapping = row.get(f"{p}_mapping", "")
    identity, coverage = _f(row.get(f"{p}_identity")), _f(row.get(f"{p}_coverage"))
    if not mapping:
        return "unmeasured (no coordinate file validated)"
    if mapping == "M4":
        return "accepted (M4)"
    if mapping == "M3":
        return "M3 (high identity and coverage, accession not confirmed)"
    if (_tf(row.get(f"{p}_accession_match")) and identity is not None and coverage is not None
            and identity >= 95 and coverage < 80):
        return "partial structure (accession confirmed, identity >= 95, coverage < 80)"
    return f"other ({mapping})"


def loss_breakdown(rows: list[dict], p: str = "method1") -> dict:
    out: dict[str, int] = {}
    for r in rows:
        if _tf(r.get(f"{p}_hit")) and r.get(f"{p}_truth_label") in POS:
            reason = loss_reason(r, p)
            out[reason] = out.get(reason, 0) + 1
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


def unmeasured_reasons(rows: list[dict], p: str = "method1") -> dict:
    out: dict[str, int] = {}
    for r in rows:
        if _tf(r.get(f"{p}_hit")) and not r.get(f"{p}_mapping"):
            reason = (r.get(f"{p}_validation_error") or "no validation_error recorded")[:90]
            out[reason] = out.get(reason, 0) + 1
    return dict(sorted(out.items(), key=lambda kv: -kv[1]))


def level_distribution(rows: list[dict]) -> dict:
    dist: dict[str, int] = {}
    for r in rows:
        level = r.get("method3_mapping") or "none"
        dist[level] = dist.get(level, 0) + 1
    return dict(sorted(dist.items()))


def analyse(rcsb_rows: list[dict], af_rows: list[dict]) -> dict:
    report: dict = {"definitions": METHODS, "rcsb": {}, "alphafold": {}}
    for method in METHODS:
        report["rcsb"][method] = {s: summarise(rcsb_rows, method, s)
                                  for s in ("all", "has_structure", "no_structure")}
    report["rcsb"]["paired_false_association_tests"] = {
        "naive1_vs_strict1": mcnemar(rcsb_rows, "method1_naive", "method1_validated_strict"),
        "naive1_vs_permissive1": mcnemar(rcsb_rows, "method1_naive", "method1_validated_permissive"),
        "naive1_vs_method3_strict": mcnemar(rcsb_rows, "method1_naive", "method3_strict"),
        "naive2_vs_strict2": mcnemar(rcsb_rows, "method2_naive", "method2_validated_strict"),
    }
    report["rcsb"]["correct_naive_hits_by_gate_outcome"] = loss_breakdown(rcsb_rows)
    report["rcsb"]["unmeasured_hit_reasons"] = unmeasured_reasons(rcsb_rows)
    report["rcsb"]["method3_level_distribution"] = level_distribution(rcsb_rows)
    report["alphafold"] = {
        "note": ("AlphaFold DB is keyed by UniProt accession, so its model is the queried entry by "
                 "construction; M4 here is close to guaranteed and is NOT evidence about validation."),
        "n": len(af_rows),
        "model_available": sum(1 for r in af_rows if r.get("method3_status") == "ok"),
        "method3_level_distribution": level_distribution(af_rows),
    }
    return report


def to_markdown(report: dict) -> str:
    lines = ["| method | stratum | queries | accepted | true | false | false-assoc. rate (95% CI) | FP/1,000 | sensitivity |",
             "|---|---|---|---|---|---|---|---|---|"]
    for method in METHODS:
        for stratum in ("all", "has_structure", "no_structure"):
            s = report["rcsb"][method][stratum]
            ci = s["false_association_rate_95ci"]
            rate = "n/a" if s["false_association_rate"] is None else f"{100*s['false_association_rate']:.1f}% ({100*ci[0]:.1f}-{100*ci[1]:.1f})"
            sens = s.get("sensitivity_correct_structure_returned")
            lines.append(f"| {method} | {stratum} | {s['n_queries']} | {s['accepted']} | {s['true_associations']} | "
                         f"{s['false_associations']} | {rate} | {s['false_associations_per_1000_queries']} | "
                         f"{'' if sens is None else f'{100*sens:.1f}%'} |")
    lines += ["", "False associations by tier (all queries):", "",
              "| method | same sequence | near-identical | different protein | unmeasured |", "|---|---|---|---|---|"]
    for method in METHODS:
        t = report["rcsb"][method]["all"]["false_association_tiers"]
        lines.append(f"| {method} | {t['same_sequence']} | {t['near_identical']} | {t['different']} | {t['unmeasured']} |")
    lines += ["", "Correct naive method-1 hits, by what the strict gate did with them:", "",
              "| outcome | n |", "|---|---|"]
    for reason, n in report["rcsb"]["correct_naive_hits_by_gate_outcome"].items():
        lines.append(f"| {reason} | {n} |")
    return "\n".join(lines)


def read(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--baselines", default=str(VALIDATION_ROOT / "results" / "baselines"))
    parser.add_argument("--out", default=str(VALIDATION_ROOT / "results" / "report"))
    args = parser.parse_args(argv)

    base = Path(args.baselines)
    rcsb_rows, af_rows = read(base / "rcsb_baseline_results.csv"), read(base / "alphafold_baseline_results.csv")
    if not rcsb_rows:
        print("No RCSB baseline results found; run 02_run_baselines.py first.", file=sys.stderr)
        return 1
    report = analyse(rcsb_rows, af_rows)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    Path(f"{out}.json").write_text(json.dumps(report, indent=2))
    Path(f"{out}.md").write_text(to_markdown(report) + "\n")
    print(to_markdown(report))
    print("\nunmeasured hit reasons:", json.dumps(report["rcsb"]["unmeasured_hit_reasons"], indent=2))
    print("\npaired tests:", json.dumps(report["rcsb"]["paired_false_association_tests"], indent=2))
    print("\nmethod3 M-level distribution (RCSB):", report["rcsb"]["method3_level_distribution"])
    print(f"\nWrote {out}.json and {out}.md", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
