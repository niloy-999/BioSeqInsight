# Manuscript plan for the BioSeqInsight 2.0 resubmission

This document maps each element of the editorial decision on
Ms. SOFTX-D-26-01242 to the specific evidence in this repository that
addresses it, and states plainly what is *not* yet evidenced.

**It is deliberately written last.** The software, the tests and the
benchmarks came first. Writing the paper first is how the v1.0 submission
ended up with claims its evidence did not support.

---

## 1. The decision being answered

> "The current work does not yet provide sufficiently strong evidence of
> software maturity, reliability, or added value over existing tools."

Three criticisms. Each needs a different kind of evidence, and conflating
them is a mistake.

| Criticism | Kind of evidence required |
|---|---|
| Maturity | Artefacts: tests, CI, packaging, documentation, versioning |
| Reliability | Measurements: failure recovery, and correctness of results |
| Added value | Comparison: against a named, documented alternative |

---

## 2. Evidence inventory

### 2.1 Maturity

| Claim | Evidence | Where |
|---|---|---|
| Layered, testable architecture | ~7 000 lines across 8 layers; `core` imports only the standard library | `src/`, `docs/architecture.md` |
| Comprehensive test suite | 432 tests, 89.6 % statement coverage (90.1 % core+structures+services) | `tests/`, `benchmarks/results/coverage_summary.json` |
| Continuous integration | 12-way OS × Python matrix, plus no-dependency, Biopython, transport, benchmark and packaging jobs | `.github/workflows/` |
| Installable | `pip install`, console scripts, optional extras, wheel installed into a clean venv in CI | `pyproject.toml`, `tests.yml:packaging` |
| Documented | README, 6 tutorials, 6 reference documents, protocols | `docs/`, `benchmarks/protocols/` |
| Maintained honestly | Every defect found during the rewrite listed with its consequence | `CHANGELOG.md`, `docs/SOFTWARE_AUDIT.md` |

**The strongest maturity argument is not the count of tests.** It is that
writing them found seven real defects in code that appeared to work,
including a reading-frame shift that produced plausible wrong peptides. Say
that, and list them. A paper that admits finding its own bugs is more
credible than one that reports a coverage percentage.

### 2.2 Reliability

Two distinct senses, and the manuscript must separate them.

**(a) Reliability under service failure.**

| Measurement | Value | Script |
|---|---|---|
| Scenarios covered | 13 | `run_fault_injection.py` |
| Requests simulated | 1 300 | |
| Requests beginning with a failure | 1 000 | |
| Recovered | 800 (**80.0 %**) | |
| Scenarios behaving as specified | 13/13 | |

Deterministic, offline, reproducible in about two seconds. This directly
answers v1.0's "9 of 25 ESM Atlas requests returned 504", which reported an
outage rather than an engineering response to one.

**(b) Reliability of the scientific output — the more important sense.**

| Claim | Evidence |
|---|---|
| Results agree with an established reference | `tests/test_biopython_parity.py`, run in a dedicated CI job |
| Values match published constants | e.g. ubiquitin 8564.74 Da vs accepted 8564.8 Da, asserted as a test |
| Retrieved structures are verified, not assumed | M0–M4 taxonomy with alignment identity, coverage and accession cross-check |
| The previous version was unreliable in exactly this way | 7 of 25 v1.0 "successes" were a different molecule |

**The 7/25 finding should be prominent in the paper, not buried.** It is the
motivation for the entire rewrite and it is a genuine contribution: it
demonstrates a failure mode that affects any tool performing unvalidated
sequence-to-structure lookup, which is most of them.

### 2.3 Added value

| Measurement | Value | Source |
|---|---|---|
| Manual workflow steps across 5 tasks | 849 | `run_workflow_comparison.py` |
| BioSeqInsight steps | 14 | |
| Reduction | 98.4 % | |
| Distinct applications | 16 → 2 | |
| Manual data transfers | 17 → 2 | |
| Capability absent from every compared tool | automatic structure identity validation | README §15 |

The baseline is published (`manual_workflow_protocol.md`) so the counts can
be disputed by someone who thinks a step is wrong.

**The honest framing of added value is not "faster".** It is: *this software
performs a verification step that the alternatives leave to the user, and
which the user usually skips.* Speed is a secondary benefit. Lead with the
verification.

---

## 3. What is NOT yet evidenced

This section exists so nothing in the manuscript outruns the repository.

| Not evidenced | Why | What the paper must say |
|---|---|---|
| Human task completion time | The usability study has not been run | Report interaction counts only; cite the protocol as future work |
| Error-rate reduction by real users | Same | Same |
| SUS or NASA-TLX scores | Same | Do not mention scores at all |
| Live retrieval success rates | Requires a dated network run | **Run `run_structure_benchmark.py` before submission and report the dated result** |
| Superiority over scripted pipelines | Out of scope by design | State the scope: users who do not write the script |

**Action required before submission:** run the live structure benchmark and
insert its output. Everything else in this repository is already measured.

---

## 4. Proposed structure (SoftwareX format)

**Title.** Keep "BioSeqInsight" but change the subtitle to foreground
validation, because that is what is new. Suggested: *"BioSeqInsight 2.0: a
desktop platform for DNA sequence analysis with validated
sequence-to-structure mapping."*

**Motivation and significance (~1 000 words).**
Open with the concrete failure: an unvalidated structure lookup returned
insulin-degrading enzyme for a preproinsulin query, and the software reported
success. Generalise: any tool that treats "a structure came back" as "the
right structure came back" has this failure mode. State the gap: existing
desktop tools do not verify identity, and the verification step is the one
users skip. Close with what this software contributes.

**Software description (~1 500 words).**
Architecture figure (the layer diagram). The M0–M4 taxonomy in a table with
its thresholds. The fallback and retry machinery. Batch processing and
projects. One paragraph on the dependency-free core and why it strengthens
rather than weakens the scientific claim.

**Illustrative examples (~800 words).**
Three: (1) the clean ubiquitin M4 case; (2) the preproinsulin case, showing
v1.0's output beside v2.0's; (3) a batch run with the mapping-level
distribution. Example 2 is the paper's centrepiece.

**Impact (~800 words).**
The three evidence categories above, with the tables. Be explicit about which
numbers are measured, which are derived, and which are not yet available.

**Limitations.** A real section, not a sentence. Reuse README §14: does not
predict structure; the propensity sketch is not a prediction; alignment is
for verification not homology search; usability results are not yet
available.

**Conclusions.** Short. Resist restating the impact section.

---

## 5. Figures and tables

| # | Content | Source |
|---|---|---|
| Fig. 1 | Layer architecture | `docs/architecture.md` |
| Fig. 2 | The M0–M4 decision path | `docs/mapping-levels.md` |
| Fig. 3 | GUI screenshot showing a mapping-level badge on a non-exact match | capture from the Structure tab |
| Fig. 4 | Fault-injection recovery by scenario (bar chart) | `benchmarks/results/fault_injection.json` |
| Fig. 5 | Batch throughput vs sequence count | `benchmarks/results/performance.json` |
| Table 1 | v1.0 vs v2.0 capability comparison | README §2 |
| Table 2 | M0–M4 definitions and thresholds | `docs/mapping-levels.md` |
| Table 3 | Live retrieval results by resource and mapping level | **run the live benchmark** |
| Table 4 | Workflow interaction counts | `benchmarks/results/workflow_comparison.json` |
| Table 5 | Comparison with existing tools | README §15 |

Figure 3 matters more than it looks: a screenshot of the software *warning*
the user that a structure is not their protein is the single clearest
statement of what this work contributes.

---

## 6. Pre-submission checklist

**Software**
- [ ] `python -m pytest -q` passes on the submission commit
- [ ] `python scripts/selfcheck.py` passes
- [ ] CI green on all matrix jobs
- [ ] Version consistent in `__init__.py`, `pyproject.toml`, `CITATION.cff`
- [ ] `CITATION.cff` author fields completed (currently placeholders)

**Evidence**
- [ ] Offline benchmarks re-run and results committed
- [ ] **Live structure benchmark run, dated, and committed**
- [ ] Every number in the manuscript traceable to a script in the repository
- [ ] No claim about human users beyond interaction counts

**Archival**
- [ ] Tagged release on GitHub
- [ ] Zenodo DOI minted and cited in the manuscript
- [ ] v1.0 preserved at `v1.0-archive`

**Manuscript**
- [ ] Limitations section present and substantive
- [ ] Structural resources cited separately (RCSB, AlphaFold, ESMFold, UniProt)
- [ ] Cover letter addresses each of the three editorial criticisms by name

---

## 7. Cover letter outline

One page. Structure:

1. Acknowledge the decision without arguing with it.
2. State that the software was rebuilt rather than revised, and say why: an
   audit of the previous version found that 7 of the 25 reported successful
   structure retrievals returned a different molecule.
3. Address each criticism in a short paragraph with a number and a pointer:
   maturity (432 tests, 89.6 % coverage, 12-way CI matrix); reliability
   (80 % recovery across 1 300 simulated failures; validated structure
   mapping); added value (849 → 14 interaction steps, and a verification step
   absent from the alternatives).
4. State what is not claimed: the usability study is specified but not run,
   and the manuscript confines itself accordingly.
5. Note that every number is reproducible from the repository, and give the
   commands.

Point 4 is the one to resist cutting. An editor who has rejected a paper for
overclaiming will read a limitations statement as a sign that the resubmission
is different in kind.
