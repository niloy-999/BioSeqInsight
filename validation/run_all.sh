#!/usr/bin/env bash
# Orchestrates PROTOCOL.md's phases in order. Safe to re-run: every step is
# cached (validation/lib/api.py on-disk cache), so interrupting and
# restarting resumes rather than re-fetching everything.
#
# Usage:
#   ./run_all.sh pilot     # Phase: ~100 proteins, sanity-check everything first
#   ./run_all.sh main      # Phase: the full ~1000-protein run (set TARGET below)
#   ./run_all.sh adversarial   # Experiment B (needs the dataset from pilot/main)
#   ./run_all.sh baselines     # Experiment C (needs the dataset)
#   ./run_all.sh threshold     # Experiment D (needs baselines + adversarial)
#   ./run_all.sh casestudies   # Experiment E part 1 (independent of the dataset)
#   ./run_all.sh report        # Experiment C statistics
#   ./run_all.sh all           # main -> adversarial -> baselines -> threshold -> casestudies -> report
#
# PROTOCOL.md Section 8: run `pilot` first and read results/pilot/*
# before committing to `main`. Do not tune the protocol based on `main`'s
# results once frozen -- only based on genuine mechanical problems found in
# `pilot`.

set -euo pipefail
cd "$(dirname "$0")"   # validation/

TARGET="${TARGET:-1000}"
PY="${PYTHON:-python3}"

DATA=data
RESULTS=results

run_pilot() {
    echo "=== PILOT: 100 proteins (sanity check before the main run) ==="
    $PY scripts/00_build_dataset.py --target 100 --out-prefix "$DATA/pilot_dataset"
    echo "Inspect $DATA/pilot_dataset_metadata.json and $DATA/pilot_dataset_accessions.csv."
    echo "Check for: empty categories, API errors above, unexpectedly slow stages."
    echo "Do not proceed to 'main' until the pilot output looks sane."
}

run_main() {
    echo "=== MAIN: $TARGET proteins (set TARGET=500 for the fallback size) ==="
    $PY scripts/00_build_dataset.py --target "$TARGET" --out-prefix "$DATA/dataset"
}

run_adversarial() {
    echo "=== Experiment B: building adversarial cases ==="
    $PY scripts/01_build_adversarial.py --dataset "$DATA/dataset_accessions.csv" \
        --out "$DATA/adversarial_cases.csv"
    echo "=== Experiment B: running adversarial cases through validation ==="
    $PY scripts/03_run_adversarial.py --cases "$DATA/adversarial_cases.csv" \
        --dataset "$DATA/dataset_accessions.csv" --out "$RESULTS/adversarial"
}

run_baselines() {
    echo "=== Experiment C: three-method baseline comparison ==="
    $PY scripts/02_run_baselines.py --dataset "$DATA/dataset_accessions.csv" \
        --out "$RESULTS/baselines"
}

run_threshold() {
    echo "=== Experiment D: identity x coverage threshold analysis ==="
    $PY scripts/04_threshold_analysis.py --baselines "$RESULTS/baselines" \
        --adversarial "$RESULTS/adversarial/adversarial_results.csv" \
        --out "$RESULTS/threshold_analysis"
}

run_casestudies() {
    echo "=== Experiment E (part 1): adjudicating the v1.0 25-protein benchmark ==="
    $PY scripts/06_case_studies.py --out "$RESULTS/case_studies"
}

run_report() {
    echo "=== Statistics: confusion matrices, CIs, McNemar (Experiment C) ==="
    $PY scripts/05_analyze_and_report.py --baselines "$RESULTS/baselines" \
        --out "$RESULTS/report"
}

case "${1:-}" in
    pilot)       run_pilot ;;
    main)        run_main ;;
    adversarial) run_adversarial ;;
    baselines)   run_baselines ;;
    threshold)   run_threshold ;;
    casestudies) run_casestudies ;;
    report)      run_report ;;
    all)
        run_main
        run_adversarial
        run_baselines
        run_threshold
        run_casestudies
        run_report
        echo "=== Done. See validation/results/ and the manuscript-figure section of README.md ==="
        ;;
    *)
        echo "Usage: $0 {pilot|main|adversarial|baselines|threshold|casestudies|report|all}" >&2
        exit 2
        ;;
esac
