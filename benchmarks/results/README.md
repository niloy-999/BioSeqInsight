# Benchmark results

Committed results from the benchmarks. Regenerate any of them with the
corresponding script in `../scripts/`.

| File | Script | Network |
|---|---|---|
| `fault_injection.json` | `run_fault_injection.py` | no |
| `performance.json` | `run_performance_benchmark.py` | no |
| `workflow_comparison.json` | `run_workflow_comparison.py` | no |
| `coverage_summary.json` | test suite with coverage | no |
| `structure_benchmark_2026-09-23_summary.json` | `run_structure_benchmark.py` | **yes** — the citable live result |

## The live structure benchmark: four runs, an audit trail

Runs are kept in order because each run after the first exists *because* the
run before it surfaced a real defect. Read them together and the history is
the evidence that these figures were checked, not asserted.

1. **`structure_benchmark_2026-09-21_prefix_summary.json`** — before any
   fix. AlphaFold retrieval: **0%** (a hard-coded, stale model-version
   guess). RCSB and ESM Atlas figures don't go through that code path and
   are representative on their own.
2. **`structure_benchmark_2026-09-22_postfix_summary_superseded.json`** —
   after the AlphaFold provider fix. AlphaFold retrieval: **98.5%**,
   confirming that fix, but `recovered_after_failure` reads 67/67 —
   every successful retrieval flagged as a recovery, which turned out to be
   a second, unrelated defect in `StructureResult.recovered` (2.0.1, item 2).
3. **`structure_benchmark_2026-09-23_summary.json`** — after the 2.0.1
   fixes. AlphaFold 98.6% (68/69), `recovered_after_failure` a plausible 0.
   Its raw output also showed ESM Atlas's `invalid_input` category (length-
   limit refusals) missing from the summary entirely, fixed and
   back-annotated onto this file (2.0.1, item 4).
4. **`structure_benchmark_2026-09-23b_summary.json`** — a second,
   independent run the same day, on 2.0.1. AlphaFold 97.1% (67/69, one more
   transient `service_error` than run 3 — ordinary day-to-day network
   variance, correctly categorised). This run's raw `best_mapping_
   distribution` summed to 68 against 69 resolved candidates: the one
   protein that got no usable structure from any provider was silently
   absent from the table instead of counted as M0 (2.0.2, item 1). Corrected
   by hand in the committed file; see its `note` for exactly what was fixed
   and what was already correct.

**This is the file to cite going forward, once produced:** a fresh run on
2.0.2 or later, where every field is native rather than hand-reconciled.
None of the four runs above needed a number changed — only added, corrected,
or annotated for a reporting gap — but a clean 2.0.2+ run is still the
right one to put in a manuscript rather than an annotated older one.

```bash
python benchmarks/scripts/run_structure_benchmark.py \
    --out benchmarks/results/structure_benchmark_$(date +%F)
```

## Checklist for any new run before citing it

Every one of these failed at some point during the sequence above — that's
why each is here:

- Every `per_provider.<name>.accounting_check` reads `"ok"`.
- `sum(best_mapping_distribution.values()) == resolved`.
- `recovered_after_failure` is small relative to `retrieved` for every
  provider — equal to it, as happened once, is the specific pattern that
  caught the `StructureResult.recovered` defect.
- `software_version` in the file matches the version whose fixes you intend
  to be citing (`CHANGELOG.md`).
- `attempted_success_rate_percent` and `success_rate_percent` are both
  present for every provider and you know which one you mean when you write
  the number down — they answer different questions (2.0.2, item 2).

## Note on the coverage figure

`coverage_summary.json` was produced by a stdlib `trace`-based measurement in
an environment without `pytest-cov`, excluding the four Tkinter view modules,
which cannot run headless. The `pytest-cov` run in CI is authoritative.

## Environment for the committed offline numbers

Performance figures depend on hardware. `performance.json` records the
Python version, platform and processor it was measured on. Numbers from a
different machine will differ in absolute terms; the *shape* of the result
(linear scaling, a per-sequence cost ratio near 1.0) should not.
