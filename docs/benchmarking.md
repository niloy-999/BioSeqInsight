# Benchmarking and reproducibility

Every number in the manuscript is produced by a script in this repository.
This document says which script, what it measures, and what it deliberately
does not measure.

---

## The four benchmarks

| Script | Network | Deterministic | Measures |
|---|---|---|---|
| `run_fault_injection.py` | no | yes | recovery from 13 scripted failure modes |
| `run_performance_benchmark.py` | no | yes (fixed seed) | throughput, scaling, alignment cost, cache latency, start-up |
| `run_workflow_comparison.py` | no | yes | interaction cost versus the documented manual workflow |
| `run_structure_benchmark.py` | **yes** | no (live services) | retrieval success and mapping levels across the curated protein set |

Three of the four need no network and give identical results on every
machine. That is the point: a reviewer can verify most of the evaluation in
under a minute without an internet connection.

---

## 1. Fault injection

```bash
python benchmarks/scripts/run_fault_injection.py
python benchmarks/scripts/run_fault_injection.py --repeats 200 --max-attempts 5
```

A scripted transport replays a defined failure pattern for every request.
Thirteen scenarios: clean success; one, two and persistent gateway timeouts;
500, 502, 503; rate limiting with `Retry-After`; 404 and 403 (which must
*not* be retried); read timeout; connection error; and a gateway answering
HTTP 200 with an HTML error page.

Reference result at 100 repeats per scenario:

| Scenario | Success | Recovered | Failed | Mean attempts |
|---|---|---|---|---|
| healthy_service | 100 | 0 | 0 | 1.00 |
| single_504 | 100 | 100 | 0 | 2.00 |
| double_504 | 100 | 100 | 0 | 3.00 |
| persistent_504 | 0 | 0 | 100 | 3.00 |
| http_500 / 502 / 503 | 100 each | 100 each | 0 | 2.00 |
| rate_limited_429 | 100 | 100 | 0 | 2.00 |
| not_found_404 | 0 | 0 | 100 | 1.00 |
| forbidden_403 | 0 | 0 | 100 | 1.00 |
| read_timeout | 100 | 100 | 0 | 2.00 |
| connection_error | 0 | 0 | 100 | 3.00 |
| html_error_page_as_200 | 100 | 100 | 0 | 2.00 |

**Aggregate: 1 000 requests began with a failure, 800 recovered — 80.0 %.**
Every scenario behaved as specified. The 200 that did not recover are the
three scenarios designed to be unrecoverable.

Backoff delays are computed by the retry policy but not slept through, so
this runs in about two seconds. Wall-clock timings are therefore not
reported here, which the script states in its output.

---

## 2. Performance

```bash
python benchmarks/scripts/run_performance_benchmark.py
python benchmarks/scripts/run_performance_benchmark.py --sizes 10,100,1000 --repeats 5
```

All inputs are generated from a fixed seed (20260921), so the workload is
identical everywhere.

Reference results (Python 3.12, x86-64 Linux):

| Measure | Value |
|---|---|
| Single 300 nt analysis | 0.54 ms |
| Single 3 000 nt analysis | 4.7 ms |
| Single 30 000 nt analysis | 50 ms |
| Single 300 aa protein | 1.6 ms |
| Batch of 10 | 761 sequences/s |
| Batch of 100 | 875 sequences/s |
| Batch of 500 | 900 sequences/s |
| Batch of 1 000 | 851 sequences/s |
| Per-sequence cost ratio, 10 → 1 000 | **1.18** (i.e. linear) |
| Cache read, 200 kB payload | 0.21 ms |
| Import + first analysis | 24.6 ms |

**Remote response times are deliberately excluded.** They are a property of
RCSB, EBI and ESM Atlas, they vary by time of day and location, and mixing
them into a scalability curve would make the curve meaningless. Live service
timings come from the structure benchmark instead.

---

## 3. Workflow comparison

```bash
python benchmarks/scripts/run_workflow_comparison.py
```

For five representative tasks, counts steps, distinct applications and manual
data transfers under the manual protocol
(`benchmarks/protocols/manual_workflow_protocol.md`) and under BioSeqInsight,
and measures the machine time BioSeqInsight takes.

| Task | Manual steps | Integrated | Reduction |
|---|---|---|---|
| T1 single-sequence characterisation | 16 | 3 | 81 % |
| T2 structure retrieval with identity check | 17 | 3 | 82 % |
| T3 batch of fifty sequences | 804 | 3 | 100 % |
| T4 reproduce an earlier analysis | 6 | 3 | 50 % |
| T5 hand the work to a collaborator | 6 | 2 | 67 % |
| **Total** | **849** | **14** | **98.4 %** |

Applications: 16 → 2. Manual data transfers: 17 → 2.

**What this is not.** Step counts are a proxy for effort, not a measurement
of human time. They are objective — the protocol is published, so two people
counting will agree — but they are not a usability result. Human completion
time and error rates require the study in
`benchmarks/protocols/usability_study_protocol.md`, which is specified but
not yet run. Nothing in this repository claims otherwise.

The manual baseline uses public web tools, because that is what the target
user uses. A bioinformatician scripting a Biopython pipeline would beat both
routes on T3, and should. The manuscript must state that scope.

---

## 4. Live structure benchmark

```bash
python benchmarks/scripts/run_structure_benchmark.py \
    --out benchmarks/results/structure_benchmark_$(date +%F)

python benchmarks/scripts/run_structure_benchmark.py --limit 5   # smoke run
```

Requires network access. Produces three files: a resolved dataset, a
per-protein-per-provider result table, and a summary.

**The dataset stores no sequence lengths.**
`benchmarks/data/benchmark_candidates.csv` holds accessions, names,
organisms, functional classes and the rationale for including each protein.
Lengths are fetched from UniProt at run time and the size strata (short
< 150 aa, medium 150–300 aa, long > 300 aa) are computed from them.

This matters: a published table that quotes a hand-typed length can silently
disagree with the sequence actually analysed. Resolving at run time makes
that impossible, and any drift in a UniProt entry between runs shows up as a
changed resolved file rather than as an invisible error.

Accessions that fail to resolve are excluded and listed by name in the
summary, so the resolved set is always auditable.

The candidate set is chosen to include hard cases on purpose: near-identical
paralogues (HBA/HBB, HRas/KRas), orthologue pairs across species, complex
partners (GroES/GroEL), intrinsically disordered proteins, membrane
receptors deposited only as fragments, and the preproinsulin case that v1.0
mis-mapped. A benchmark composed only of ubiquitin and lysozyme would report
a flattering number and test nothing.

---

## Reproducing the published figures

```bash
git clone https://github.com/niloy-999/BioSeqInsight.git
cd BioSeqInsight
pip install -e ".[validation]"

python -m pytest -q --cov=bioseqinsight            # 432 tests, coverage
python benchmarks/scripts/run_fault_injection.py    # reliability
python benchmarks/scripts/run_performance_benchmark.py
python benchmarks/scripts/run_workflow_comparison.py
python benchmarks/scripts/run_structure_benchmark.py --out results/live   # needs network
```

The first four take about a minute in total. The last depends on how fast the
public services respond.

Results land in `benchmarks/results/` as JSON and CSV. The committed files
are the published ones; dated runs are gitignored so a local run cannot
quietly overwrite the record.

---

## Honest reporting checklist

Before any number goes into the manuscript:

- [ ] Which script produced it, and is that script in the repository?
- [ ] Does it need network access? If so, is the run date recorded?
- [ ] Is it deterministic? If not, how many repeats, and is variance reported?
- [ ] Is it a measurement or a derivation? Step counts are derived; throughput
      is measured; human task time is neither until the study is run.
- [ ] Does the claim match what was measured, or is it broader?
- [ ] Are the failure cases reported alongside the successes?

The last two are the ones v1.0 failed. "25/25 structures retrieved" was a
count of HTTP 200 responses presented as a count of correct answers.
