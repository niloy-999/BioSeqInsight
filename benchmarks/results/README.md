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

## The live structure benchmark: five runs, an audit trail, now closed

Runs are kept in order because each run after the first exists *because* the
run before it surfaced a real defect — until the fifth, which surfaced none.

1. **`structure_benchmark_2026-09-21_prefix_summary.json`** — before any
   fix. AlphaFold retrieval: **0%** (a hard-coded, stale model-version
   guess). RCSB and ESM Atlas figures don't go through that code path and
   are representative on their own.
2. **`structure_benchmark_2026-09-22_postfix_summary_superseded.json`** —
   after the AlphaFold provider fix. AlphaFold retrieval: **98.5%**,
   confirming that fix, but `recovered_after_failure` reads 67/67 —
   a second, unrelated defect in `StructureResult.recovered` (2.0.1, item 2).
3. **`structure_benchmark_2026-09-23_summary.json`** — after the 2.0.1
   fixes. AlphaFold 98.6% (68/69). Raw output was missing the
   `invalid_input` category for ESM Atlas, back-annotated (2.0.1, item 4).
4. **`structure_benchmark_2026-09-23b_summary.json`** — a second run the
   same day, still on 2.0.1. AlphaFold 97.1% (67/69, one more transient
   `service_error` — ordinary variance). `best_mapping_distribution` summed
   to 68 against 69 resolved: one genuinely-M0 protein was silently absent
   from the table instead of counted, corrected by hand (2.0.2, item 1).
5. **`structure_benchmark_2026-09-23c_summary.json` — THE CITABLE RESULT.**
   A third run the same day, on 2.0.2. AlphaFold 98.6% (68/69), overall
   exact-mapping rate 98.6%, any-structure rate 100%. **Every field is
   native — nothing added, corrected, or back-annotated.** Every
   `accounting_check` reads `"ok"` on its own; `best_mapping_distribution`
   sums to `resolved` on its own; strata `exact_matches` (16+22+30=68)
   independently cross-validates against `best_mapping_distribution.M4`
   (68). This is the first run in the sequence that needed no correction at
   all, which is what makes it the one to put in the manuscript.

**The audit trail is closed as of run 5.** A future re-run is still worth
doing closer to submission for a fresh timestamp, but the defect-finding
part of this process — the reason four earlier runs exist — is done: two
real provider/model defects (items 1-2, 2.0.1) and two real
benchmark-reporting defects (item 4, 2.0.1; item 1, 2.0.2) were found this
way, none of them reachable from the offline test suite alone, all of them
now covered by regression tests so they cannot recur silently.

```bash
python benchmarks/scripts/run_structure_benchmark.py \
    --out benchmarks/results/structure_benchmark_$(date +%F)
```

## Checklist for any future run before citing it

Every one of these failed at some point during runs 1-4 above — that's why
each is here, and why run 5 passing all of them is what ends the audit:

- Every `per_provider.<name>.accounting_check` reads `"ok"`.
- `sum(best_mapping_distribution.values()) == resolved`.
- Strata `exact_matches` sums to `best_mapping_distribution["M4"]`.
- `recovered_after_failure` is small relative to `retrieved` for every
  provider — equal to it, as happened once, is the specific pattern that
  caught the `StructureResult.recovered` defect.
- `software_version` in the file matches the version whose fixes you intend
  to be citing (`CHANGELOG.md`).
- `attempted_success_rate_percent` and `success_rate_percent` are both
  present and you know which one you mean when you write the number down —
  they answer different questions (2.0.2, item 2).

## Note on the coverage figure

`coverage_summary.json` was produced by a stdlib `trace`-based measurement in
an environment without `pytest-cov`, excluding the four Tkinter view modules,
which cannot run headless. The `pytest-cov` run in CI is authoritative.

## Environment for the committed offline numbers

Performance figures depend on hardware. `performance.json` records the
Python version, platform and processor it was measured on. Numbers from a
different machine will differ in absolute terms; the *shape* of the result
(linear scaling, a per-sequence cost ratio near 1.0) should not.
