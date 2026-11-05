# BioSeqInsight validation package

Implements the six-experiment plan (A-F) for the Bioinformatics Advances
submission: a stratified ~1,000-protein benchmark, an adversarial benchmark,
a three-method baseline comparison, an identity x coverage threshold
analysis, real-data case studies (including adjudication of the v1.0
25-protein benchmark), and statistics with confidence intervals. Experiment
F (software engineering evidence: tests, fault injection, performance,
workflow comparison) is already complete in `../benchmarks/`, not here.

**Everything in this package is new code, written and offline-tested in
this session but never run against the live internet** (the environment
that wrote it has no network access). Every pure-logic component --
dataset classification, stratified sampling, adversarial case construction,
ground-truth adjudication, PDB-ID resolution, scoring, and every statistic
in `lib/stats.py` -- has been verified against synthetic/mocked data with
assertions that check the actual numbers, not just "it runs." What has
**not** been exercised is the real network calls themselves: real UniProt
query syntax, real RCSB/SIFTS response shapes, real rate limits. Run the
pilot (Section 2 below) first, specifically to catch anything that only
shows up against the real APIs.

## 1. Why UniProt (not NCBI/RefSeq), and InterPro only as a refinement

Every resource `bioseqinsight` queries -- RCSB's own cross-reference check,
AlphaFold DB, the canonical-sequence fetch itself -- is keyed by **UniProt
accession**. Building the benchmark from NCBI/RefSeq accessions would add an
error-prone identifier-mapping step before the software under test could
even be invoked, and would not match what a real user of this software
actually queries with. So:

- **UniProt REST API** is the dataset's primary source (`lib/api.py:
  uniprot_search`, `uniprot_entry`).
- **RCSB PDB** (`rcsb_sequence_search`, `rcsb_polymer_entity`,
  `rcsb_entry`, raw file download) supplies candidate structures and one of
  the two independent ground-truth signals.
- **PDBe SIFTS** (`sifts_uniprot_mappings`) supplies the other independent
  ground-truth signal and the ``no manual curation'' residue-level mapping
  used to classify fragments/partial structures in Experiment A.
- **InterPro is not used.** Pfam family grouping (fetched directly from each
  UniProt entry's own cross-references) is what `01_build_adversarial.py`
  uses to find paralogue pairs; this needed only one extra field already
  present in the UniProt entry JSON, so a separate InterPro integration
  was not necessary. If you want InterPro-based family grouping instead
  (e.g. for finer-grained domain-level paralogue detection), the
  Pfam-grouping block in `01_build_adversarial.py` (`fetch_pfam_ids`) is the
  one function to swap.

Nothing here requires you to hand-pick or hand-download a single accession.
`00_build_dataset.py` queries, classifies, and samples automatically;
`01_build_adversarial.py` constructs adversarial pairs automatically from
the result.

## 2. Running it: phase order

```bash
cd validation
pip install -e ".[http]" --config-settings editable_mode=compat 2>/dev/null || true  # from the repo root, if not already installed
./run_all.sh pilot          # ~100 proteins; read the output before continuing
# inspect data/pilot_dataset_metadata.json: category counts, errors, timing
./run_all.sh all            # TARGET=1000 by default; TARGET=500 ./run_all.sh all for the fallback size
```

Or step by step (`./run_all.sh <phase>` for any individual phase — see the
script header for the full list). Every script also runs standalone with
`--help` for its exact arguments.

**Expected runtime**, extrapolated from `benchmarks/scripts/
run_structure_benchmark.py`'s measured ~8-9 s/protein against three live
structural resources (`../benchmarks/results/structure_benchmark_2026-09-23_summary.json`):
a 1,000-protein main run, querying UniProt + RCSB + SIFTS per candidate plus
three retrieval methods in Experiment C, is plausibly **several hours**, not
minutes. Budget accordingly, and definitely run the 100-protein pilot first
to get a real per-protein timing figure on your own connection before
committing to the full run.

## 3. What each script produces, and where it feeds the manuscript

| Script | Experiment | Key output | Manuscript placeholder it fills |
|---|---|---|---|
| `00_build_dataset.py` | A | `data/dataset_accessions.csv`, `data/dataset_metadata.json` | Table "benchmark dataset composition"; Figure 3 in the original plan |
| `01_build_adversarial.py` | B (construction) | `data/adversarial_cases.csv` | Table "adversarial benchmark composition" |
| `03_run_adversarial.py` | B (execution) | `results/adversarial/adversarial_results.csv`, `adversarial_summary.json` | Figure "adversarial validation" / confusion matrix |
| `02_run_baselines.py` | C | `results/baselines/{rcsb,alphafold}_baseline_results.csv` | raw data for the baseline comparison table |
| `05_analyze_and_report.py` | C (statistics) | `results/report.json`: confusion matrices, precision/recall/F1/FAR with Wilson CIs, McNemar p-values, per resource | Table "baseline comparison" (the plan's Figure 4) |
| `04_threshold_analysis.py` | D | `results/threshold_analysis.json`: AUC(identity) vs AUC(identity x coverage), 2D grid | Figure "identity x coverage discrimination" (the plan's Section 8 figure) |
| `06_case_studies.py` | E (part 1) | `results/case_studies.csv`, `case_studies_summary.json` | Table~discordance and the case-study prose in `manuscript/sections/05_evaluation.tex` (search for `PENDING`) |

Two real-data worked-example transcripts for Experiment E (one fragment
case, one clean match) are not scripted separately -- once `case_studies.csv`
exists, pick one `fragment` row and one `exact_full_length` row with
`mapping_2026 == "M3"` / `"M4"` respectively and run
`bioseqinsight structure <accession> --all` by hand for the transcript; this
is the one place in the pipeline a human picks a representative example
rather than everything being scripted, and it should stay that way so a
reader can tell the example wasn't cherry-picked by an automated search for
the most favourable case.

## 4. The corrected v1.0 case-study framing

The 7 "accession-discordant" queries in v1.0's own Table 3 are **not**
confirmed wrong-molecule results -- they are flagged discordances v1.0 could
not resolve, because it had no sequence comparison step at all, and because
its reported confidence values are **AlphaFold pLDDT** (confidence in the
model AlphaFold actually returned), not identity to the query. Do not
describe these as "7 errors" anywhere you write about this; describe them as
"7 discordant cases, adjudicated as follows" and cite
`results/case_studies.csv`. See `../docs/SOFTWARE_AUDIT.md` Section 2 for
the full, precise account and the exact Table 3 values quoted verbatim.

## 5. Protocol integrity

`PROTOCOL.md` is frozen. If running the pilot surfaces a genuine
methodological problem (not just a disappointing number), fix it, log the
change in `PROTOCOL.md`'s Amendments section with a reason, and re-run --
do not silently edit the criteria and re-run without a logged reason.

## 6. Known gaps to be aware of before you run this

- `uniprot_query_for_stratum`'s query string
  (`reviewed:true AND length:[X TO Y]`) is standard, documented UniProtKB
  REST syntax, but has not been executed against the live endpoint. If it
  returns zero results or an error, check
  <https://www.uniprot.org/help/query-fields> and fix that one line; every
  downstream step classifies from the returned JSON, not from further
  query-string assumptions.
- `rcsb_sequence_search`'s request body follows RCSB's documented search API
  schema (`https://search.rcsb.org/rcsbsearch/v2/query`) but has likewise
  not been exercised live. If RCSB's schema has changed since this was
  written, the error will surface immediately on the first pilot query
  rather than silently producing wrong results, because `FetchError` is
  raised (not swallowed) on a non-2xx response other than 404.
- Method 1 vs Method 2 for RCSB (PROTOCOL.md Section 4) differ only in
  whether `identity_cutoff` is passed to RCSB's own search request --
  both use RCSB's own ranking unfiltered otherwise. If you want a
  genuinely different "naive" baseline (e.g. literally the first hit with
  no cutoff at all, including very low-identity hits RCSB itself would
  normally filter), set `evalue_cutoff` higher in `method1_rcsb`.
- B3 (wrong isoform) cases depend on genes in your sampled dataset actually
  having a solved structure for a *non-canonical* isoform specifically --
  this may be a small number even in a 1,000-protein sample. The pilot's
  `data/pilot_dataset_metadata.json` will show how many isoform-bearing
  candidates were drawn; if too few produce resolvable B3 cases in the
  main run, that is worth noting as a limitation rather than forcing the
  count artificially.
