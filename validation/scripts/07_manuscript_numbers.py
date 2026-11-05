#!/usr/bin/env python3
# AMENDMENT-9
"""Compute every validation-study number the manuscript quotes.

The manuscript contains NO hand-typed validation-study numbers: each is looked
up from ``numbers.tex`` (macro ``\\N{key}``) which this script writes from the
results files, so the paper cannot disagree with the data. ``numbers.json``
carries the same values plus the raw tables the figures are drawn from.

Reads (all produced by the earlier scripts):
  results/report.json                  05_analyze_and_report.py
  results/threshold_analysis.json      04_threshold_analysis.py
  results/adversarial/adversarial_results.csv   03_run_adversarial.py
  results/case_studies_summary.json    06_case_studies.py
  data/dataset_metadata.json           00_build_dataset.py

Writes results/manuscript/numbers.json and numbers.tex.

Usage
-----
    python validation/scripts/07_manuscript_numbers.py
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

METHODS = {
    "n1": "method1_naive", "s1": "method1_validated_strict", "p1": "method1_validated_permissive",
    "n2": "method2_naive", "s2": "method2_validated_strict", "p2": "method2_validated_permissive",
    "s3": "method3_strict", "p3": "method3_permissive",
}
STRATA = {"all": "all", "has": "has_structure", "no": "no_structure"}
CLASSES = ("B1", "B2", "B3", "B4", "B5")
LEVELS = ("M0", "M1", "M2", "M3", "M4")


def pct(value, digits: int = 1) -> str:
    return "n/a" if value is None else f"{100 * value:.{digits}f}"


def num(value, digits: int = 1) -> str:
    return "n/a" if value is None else f"{value:.{digits}f}"


def pval(p) -> str:
    """LaTeX-safe p-value."""
    if p is None:
        return "n/a"
    if p >= 0.001:
        return f"{p:.3f}"
    mantissa, exponent = f"{p:.1e}".split("e")
    return f"{mantissa}\\times10^{{{int(exponent)}}}"


def class_of(label: str) -> str:
    return label.split("_", 1)[0]


def adversarial_summary(rows: list[dict]) -> dict:
    out: dict = {"total": len(rows), "scored": 0, "classes": {}}
    by_class: dict[str, list[dict]] = collections.defaultdict(list)
    for r in rows:
        by_class[class_of(r.get("class", ""))].append(r)
        out["scored"] += r.get("status") == "ok"
    out["unresolved"] = out["total"] - out["scored"]
    for cls in CLASSES:
        group = by_class.get(cls, [])
        ok = [r for r in group if r.get("status") == "ok"]
        levels = collections.Counter(r.get("mapping_level") for r in ok)
        block = {
            "built": len(group), "n": len(ok),
            "levels": {lv: levels.get(lv, 0) for lv in LEVELS},
            "partial": sum(1 for r in ok if str(r.get("partial_structure")).lower() == "true"),
            "accession_confirmed": sum(1 for r in ok if str(r.get("accession_match")).lower() == "true"),
            "m4_cases": [{"query": r.get("query_accession"), "pdb": r.get("resolved_pdb_id"),
                          "identity": r.get("identity_percent"), "coverage": r.get("coverage_percent"),
                          "accession_match": r.get("accession_match")} for r in ok if r.get("mapping_level") == "M4"],
        }
        if cls == "B3":
            valid = []
            for r in ok:
                try:
                    if float(r["forced_structure_isoform_length"]) < 0.95 * float(r["query_canonical_length"]):
                        valid.append(r)
                except (KeyError, TypeError, ValueError):
                    pass
            block["valid_rule_n"] = len(valid)
            block["valid_rule_M4"] = sum(1 for r in valid if r.get("mapping_level") == "M4")
            block["valid_rule_available"] = bool(ok) and all(
                r.get("query_canonical_length") not in (None, "") for r in ok)
        out["classes"][cls] = block
    return out


def adversarial_diagnostics(rows: list[dict]) -> dict:
    """Facts established by inspecting the frozen adversarial results (Amendment 9)."""
    out: dict[str, str] = {}
    ok = [r for r in rows if r.get("status") == "ok"]

    def f(r, key):
        try:
            return float(r[key])
        except (KeyError, TypeError, ValueError):
            return None

    # B3: cases whose isoform is >5% shorter than the canonical sequence and that reached M4.
    iso, sub = [], []
    for r in ok:
        if class_of(r["class"]) != "B3" or r.get("mapping_level") != "M4":
            continue
        i, c, s_, q = (f(r, "forced_structure_isoform_length"), f(r, "query_canonical_length"),
                       f(r, "subject_length"), f(r, "query_length"))
        if None in (i, c, s_, q) or not i < 0.95 * c:
            continue
        iso.append(100 * i / c)
        sub.append(100 * s_ / q)
    for key, vals, fn in (("iso-min", iso, min), ("iso-max", iso, max), ("sub-min", sub, min), ("sub-max", sub, max)):
        out[f"adv-B3-valid-M4-{key}"] = num(fn(vals)) if vals else "n/a"

    # B2: confirmed accession but not flagged partial.
    unflagged = [r for r in ok if class_of(r["class"]) == "B2" and str(r.get("accession_match")).lower() == "true"
                 and str(r.get("partial_structure")).lower() != "true"]
    ids = [f(r, "identity_percent") for r in unflagged if f(r, "identity_percent") is not None]
    out["adv-B2-unflag-n"] = str(len(unflagged))
    out["adv-B2-unflag-lowid"] = str(sum(1 for v in ids if v < 95))
    out["adv-B2-unflag-highcov"] = str(sum(1 for r in unflagged if (f(r, "coverage_percent") or 0) >= 80))
    out["adv-B2-unflag-id-min"] = num(min(ids)) if ids else "n/a"
    out["adv-B2-unflag-id-max"] = num(max(ids)) if ids else "n/a"

    # B1: the pair(s) that reached M4.
    b1 = [r for r in ok if class_of(r["class"]) == "B1" and r.get("mapping_level") == "M4"]
    for key, field in (("sub", "subject_length"), ("q", "query_length")):
        v = f(b1[0], field) if b1 else None
        out[f"adv-B1-m4-{key}"] = str(int(v)) if v is not None else "n/a"
    return out


def compute(report: dict, threshold: dict, adversarial_rows: list[dict], cases: dict,
            dataset: dict) -> tuple[dict, dict]:
    flat: dict[str, str] = {}
    rc = report["rcsb"]

    for short, method in METHODS.items():
        for st_short, stratum in STRATA.items():
            s = rc[method][stratum]
            p = f"{short}-{st_short}"
            flat[f"{p}-q"] = str(s["n_queries"])
            flat[f"{p}-acc"] = str(s["accepted"])
            flat[f"{p}-true"] = str(s["true_associations"])
            flat[f"{p}-false"] = str(s["false_associations"])
            flat[f"{p}-per1000"] = num(s["false_associations_per_1000_queries"])
            flat[f"{p}-far"] = pct(s["false_association_rate"])
            ci = s["false_association_rate_95ci"]
            flat[f"{p}-far-lo"] = pct(ci[0]) if ci else "n/a"
            flat[f"{p}-far-hi"] = pct(ci[1]) if ci else "n/a"
            if "sensitivity_correct_structure_returned" in s:
                flat[f"{p}-sens"] = pct(s["sensitivity_correct_structure_returned"])
                flat[f"{p}-sens-lo"] = pct(s["sensitivity_95ci"][0])
                flat[f"{p}-sens-hi"] = pct(s["sensitivity_95ci"][1])
            tiers = s.get("false_association_tiers") or {}
            for t_short, t_key in (("same", "same_sequence"), ("near", "near_identical"),
                                   ("diff", "different"), ("unm", "unmeasured")):
                flat[f"{p}-tier-{t_short}"] = str(tiers[t_key]) if t_key in tiers and tiers[t_key] is not None else "n/a"

    gate = {"accepted": 0, "partial": 0, "unmeasured": 0, "m3": 0, "other": 0}
    for reason, n in rc["correct_naive_hits_by_gate_outcome"].items():
        key = ("accepted" if reason.startswith("accepted") else "partial" if reason.startswith("partial")
               else "unmeasured" if reason.startswith("unmeasured") else "m3" if reason.startswith("M3") else "other")
        gate[key] += n
    for k, v in gate.items():
        flat[f"gate-{k}"] = str(v)
    flat["gate-total"] = str(sum(gate.values()))
    flat["gate-rejected"] = str(sum(gate.values()) - gate["accepted"])

    tests = rc["paired_false_association_tests"]
    for key, name in (("naive1_vs_strict1", "n1-s1"), ("naive1_vs_permissive1", "n1-p1"),
                      ("naive1_vs_method3_strict", "n1-s3"), ("naive2_vs_strict2", "n2-s2")):
        flat[f"p-{name}"] = pval(tests[key]["p_value"])
        flat[f"disc-{name}"] = str(tests[key]["false_assoc_only_in_first"])

    for lv in LEVELS:
        flat[f"lvl-rcsb-{lv}"] = str(rc["method3_level_distribution"].get(lv, 0))
        flat[f"lvl-af-{lv}"] = str(report["alphafold"]["method3_level_distribution"].get(lv, 0))
    flat["af-n"] = str(report["alphafold"]["n"])
    flat["af-avail"] = str(report["alphafold"]["model_available"])

    for k, v in threshold["thresholds_used"].items():
        flat[f"th-{k.replace('_', '-')}"] = f"{v:g}"
    for tag, key in (("n", "naive_positive_vs_naive_negative"), ("a", "naive_positive_vs_adversarial_negative")):
        a = threshold[key]
        flat[f"auc-{tag}-npos"], flat[f"auc-{tag}-nneg"] = str(a["n_positive"]), str(a["n_negative"])
        flat[f"auc-{tag}-id"] = f"{a['auc_identity']:.3f}"
        flat[f"auc-{tag}-cov"] = f"{a['auc_coverage']:.3f}"
        flat[f"auc-{tag}-ic"] = f"{a['auc_identity_x_coverage']:.3f}"
    ab = threshold["ablation_at_M4_thresholds"]
    for src, dst in (("negatives_total", "neg"), ("pass_identity_and_coverage_only", "neg-ic"),
                     ("pass_identity_coverage_and_accession", "neg-ica"), ("positives_total", "pos"),
                     ("positives_pass_identity_and_coverage", "pos-ic"), ("positives_pass_all_three", "pos-all")):
        flat[f"abl-{dst}"] = str(ab[src])
    nf = threshold["unrelated_pair_noise_floor"]
    flat["nf-n"] = str(nf["n_unrelated_pairs"])
    for q in ("p50", "p90", "p95", "p99"):
        flat[f"nf-{q}-all"] = num(nf["identity_percentiles_all_pairs"][q])
        cond = nf["identity_percentiles_given_coverage_ge_related"]
        flat[f"nf-{q}-cov"] = num(cond[q]) if cond else "n/a"
    flat["nf-n-cov"] = str(nf["n_with_coverage_at_least_related_threshold"])
    flat["nf-hilow"] = str(nf["high_identity_low_coverage_pairs"])
    flat["nf-m2"] = str(nf["reach_M2_criteria"])
    flat["nf-m2-pct"] = pct(nf["reach_M2_criteria_fraction"])
    flat["nf-id50"] = str(nf["reach_identity_50_and_related_coverage"])
    flat["nf-high"] = str(nf["reach_high_criteria"])
    flat["nf-exact"] = str(nf["reach_exact_identity_and_coverage"])
    fr = threshold["fragments_B2"]
    flat["frag-n"], flat["frag-conf"] = str(fr["n"]), str(fr.get("accession_confirmed", "n/a"))
    flat["frag-med-id"], flat["frag-med-cov"] = num(fr["median_identity"]), num(fr["median_coverage"])

    adv = adversarial_summary(adversarial_rows)
    flat["adv-total"], flat["adv-scored"], flat["adv-unresolved"] = str(adv["total"]), str(adv["scored"]), str(adv["unresolved"])
    for cls, block in adv["classes"].items():
        flat[f"adv-{cls}-built"], flat[f"adv-{cls}-n"] = str(block["built"]), str(block["n"])
        flat[f"adv-{cls}-partial"], flat[f"adv-{cls}-accconf"] = str(block["partial"]), str(block["accession_confirmed"])
        for lv in LEVELS:
            flat[f"adv-{cls}-{lv}"] = str(block["levels"][lv])
    b4 = adv["classes"]["B4"]
    flat["adv-B4-M2-pct"] = pct(b4["levels"]["M2"] / b4["n"]) if b4["n"] else "n/a"
    b3 = adv["classes"]["B3"]
    avail = b3.get("valid_rule_available")
    flat["adv-B3-valid-n"] = str(b3["valid_rule_n"]) if avail else "n/a"
    flat["adv-B3-valid-M4"] = str(b3["valid_rule_M4"]) if avail else "n/a"
    b1m4 = adv["classes"]["B1"]["m4_cases"]
    for key, field in (("id", "identity"), ("cov", "coverage")):
        try:
            flat[f"adv-B1-m4-{key}"] = num(float(b1m4[0][field])) if b1m4 else "n/a"
        except (TypeError, ValueError):
            flat[f"adv-B1-m4-{key}"] = "n/a"
    flat["adv-nonB3-M4"] = str(sum(adv["classes"][c]["levels"]["M4"] for c in ("B1", "B2", "B4", "B5")))
    flat["adv-nonB3-n"] = str(sum(adv["classes"][c]["n"] for c in ("B1", "B2", "B4", "B5")))

    flat.update(adversarial_diagnostics(adversarial_rows))

    for tag, key in (("disc", "v1_discordant_queries"), ("oth", "v1_other_queries")):
        c = cases[key]
        for src, dst in (("n", "n"), ("top_hit_found", "found"), ("top_hit_entry_contains_query_accession", "contains"),
                         ("top_hit_entry_does_not_contain_query_accession", "notcontains"),
                         ("top_hit_ambiguous_or_unknown", "amb"), ("top_hit_validated_M4", "M4"),
                         ("top_hit_validated_M3", "M3"), ("top_hit_validated_below_M3", "low")):
            flat[f"cs-{tag}-{dst}"] = str(c[src])

    flat["ds-n"] = str(dataset["achieved"])
    flat["ds-orgs"] = str(dataset["sampling_diagnostics"]["distinct_organisms"])
    flat["ds-toporg-pct"] = pct(dataset["sampling_diagnostics"]["top_organism_fraction"])
    notes = " ".join(dataset["sampling_diagnostics"].get("sampling_notes", []))
    m = re.search(r"release (\d{4}_\d{2})", notes)
    flat["ds-release"] = m.group(1).replace("_", "\\_") if m else "n/a"
    for cat, short in (("exact_full_length", "exact"), ("fragment", "fragment"), ("isoform_bearing", "isoform"),
                       ("paralogue_family", "paralogue"), ("no_structure", "nostruct")):
        flat[f"ds-cat-{short}"] = str(dataset["category_counts"].get(cat, "n/a"))
    flat["ds-has"], flat["ds-no"] = flat["n1-has-q"], flat["n1-no-q"]

    raw = {
        "report_rcsb": {m_: {s: {k: rc[m_][s][k] for k in ("accepted", "true_associations", "false_associations", "false_association_tiers")}
                             for s in STRATA.values()} for m_ in METHODS.values()},
        "gate_outcomes": rc["correct_naive_hits_by_gate_outcome"],
        "adversarial": adv,
        "grid_positive": threshold.get("grid_positive"), "grid_negative": threshold.get("grid_negative"),
        "identity_bins": threshold.get("identity_bins"), "coverage_bins": threshold.get("coverage_bins"),
        "thresholds": threshold["thresholds_used"],
    }
    return flat, raw


def to_tex(flat: dict) -> str:
    lines = ["% Generated by validation/scripts/07_manuscript_numbers.py -- do not edit by hand.",
             "\\providecommand{\\N}[1]{\\csname num:#1\\endcsname}"]
    for key in sorted(flat):
        lines.append(f"\\expandafter\\gdef\\csname num:{key}\\endcsname{{{flat[key]}}}")
    return "\n".join(lines) + "\n"


def load_json(path: Path) -> dict:
    if not path.is_file():
        sys.exit(f"missing {path} -- run the earlier phases first")
    return json.loads(path.read_text(encoding="utf-8"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--results", default=str(ROOT / "results"))
    parser.add_argument("--dataset-metadata", default=str(ROOT / "data" / "dataset_metadata.json"))
    parser.add_argument("--out", default=str(ROOT / "results" / "manuscript"))
    args = parser.parse_args(argv)
    res = Path(args.results)

    report = load_json(res / "report.json")
    threshold = load_json(res / "threshold_analysis.json")
    cases = load_json(res / "case_studies_summary.json")
    dataset = load_json(Path(args.dataset_metadata))
    adv_path = res / "adversarial" / "adversarial_results.csv"
    if not adv_path.is_file():
        sys.exit(f"missing {adv_path}")
    with open(adv_path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    flat, raw = compute(report, threshold, rows, cases, dataset)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "numbers.json").write_text(json.dumps({"flat": flat, "raw": raw}, indent=2))
    (out / "numbers.tex").write_text(to_tex(flat))
    print(f"{len(flat)} numbers written to {out}/numbers.json and numbers.tex")
    for key in ("ds-n", "n1-all-false", "s1-all-false", "s1-has-sens", "adv-total", "adv-scored",
                "adv-B1-n", "adv-B1-M4", "adv-B2-partial", "adv-B2-n", "adv-B3-valid-n", "adv-B3-valid-M4",
                "nf-m2", "gate-partial", "gate-unmeasured"):
        print(f"  {key:18s} {flat.get(key)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
