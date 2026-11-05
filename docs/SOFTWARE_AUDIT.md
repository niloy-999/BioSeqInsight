# Audit of BioSeqInsight 1.0

This document records what a systematic review of version 1.0 found. It is
published because the central design decision in version 2.0 — automatic
validation of sequence-to-structure mapping — only makes sense once the
problem it solves is visible.

It is also here because the editorial decision on the v1.0 manuscript said
the work did not "provide sufficiently strong evidence of software maturity,
reliability, or added value over existing tools." The most useful response to
that is not a rebuttal but an inventory.

---

## 1. Scope

The audit covered the v1.0 source (1 242 lines across three modules), its
test file (108 lines), its README, and the benchmark table reported in the
submitted manuscript (25 proteins, RCSB / AlphaFold DB / ESM Atlas).

---

## 2. The finding that mattered

The v1.0 manuscript's Table 3 ("Individual results for the 25-protein
benchmark evaluated on 28 August 2026") reported, for every query: sequence
length, whether the live ESM Atlas request returned a model or an HTTP 504,
the mean per-residue pLDDT of that ESM Atlas model when it succeeded, the
pLDDT of the retrieved AlphaFold DB model, and a note. The note column
flagged 7 of the 25 queries with **"mapped \<accession\>"**, meaning the
retrieval step had associated the query with a UniProt accession other than
the one queried, reproduced verbatim below:

| Query | Query length (aa) | AF pLDDT | Note (v1.0's own flag) |
|---|---|---|---|
| P00974 (BPTI) | 100 | 97.80 | mapped Q52V24 |
| P0A6F9 (GroES) | 97 | 90.80 | mapped GroEL |
| P61769 ($\beta_2$-microglobulin) | 119 | 86.82 | mapped Q30597 |
| P01308 (insulin) | 110 | 93.99 | mapped insulin-degrading enzyme (IDE) |
| P68871 (Hb-$\beta$) | 147 | 98.06 | mapped Hb-$\alpha$ |
| P04637 (p53) | 393 | 75.99 | mapped RPB1 |
| UBIQ_76 (ubiquitin, 76 aa) | 76 | 88.65 | mapped P0CG48 |

**A precise, important point about what these numbers are.** The "AF pLDDT"
column is AlphaFold's own per-residue confidence in *its own predicted
structure* for whatever accession the retrieval step landed on. It is not a
sequence-identity or coverage measurement between the query and the
returned structure, and v1.0's pipeline did not compute one — nothing in
Table 3 says how similar P00974's sequence is to Q52V24's. A high pLDDT
(90-98% for five of the seven flagged cases above) meant AlphaFold was
confident in the structure it returned; it said nothing about whether that
structure was for the queried protein. That is precisely the gap v2.0's
identity/coverage/accession check (this document's Section 4;
manuscript Section on the M0-M4 classification) closes, and precisely why a
high-confidence prediction for the wrong protein is one of the failure
modes the M0-M4 scale is designed to catch (a case that would present as a
strong AlphaFold pLDDT with a weak or absent sequence match is exactly the
"M1, high pLDDT but low identity" pattern).

The seven flagged cases are not, on this evidence alone, seven confirmed
wrong molecules — v1.0's own manuscript did not claim that; it explicitly
stated a returned hit "does not by itself establish that the retrieved
structure is the exact biological target." What v1.0 lacked was any
mechanism to go further. The plausible outcomes differ case by case: the
ubiquitin case plausibly maps to P0CG48 because P0CG48 (polyubiquitin-C)
literally contains the 76-residue ubiquitin sequence as a repeat unit; the
GroES/GroEL and Hb-$\alpha$/Hb-$\beta$ cases are plausibly multi-chain
entries containing both partners; the insulin/IDE, $\beta_2$m, BPTI, and
p53 cases have no such obvious benign explanation from length and name
alone. **These are hypotheses, not findings.** They are adjudicated
properly, with live sequence alignment against each retrieved chain and
independent PDBe SIFTS / RCSB cross-reference agreement (not
BioSeqInsight's own alignment used as its own ground truth), by
`validation/scripts/06_case_studies.py`. Its output
(`validation/results/case_studies.csv`) is the authoritative record; this
table is the historical primary source it starts from.

**What v2.0 does about the underlying gap.** Every retrieval is aligned
against the query and cross-checked against the accession, and the result
is reported on the M0-M4 scale, so a pLDDT number is never presented
without the identity and coverage measurements needed to interpret it.

---

## 3. Architectural findings

| Finding | Consequence | Status in 2.0 |
|---|---|---|
| Every analysis function returned a preformatted display string | batch processing, export and reproducibility were all impossible without re-parsing text | typed dataclasses with `.to_dict()` |
| Bare-module imports (`import sequence_operations`) | the program only ran from its own directory | package-relative imports, installable |
| GUI and computation interleaved | the scientific code could not be tested without a display server | `core` has no GUI or I/O imports |
| No packaging metadata | no `pip install`, no pinned dependencies, no entry points | `pyproject.toml`, console scripts, extras |
| No CI | nothing was verified on any platform but the author's | 12-way matrix plus lint, packaging and benchmark jobs |
| Retry logic widened the timeout and tried three times | a 504 at second 1 waited out the full extended timeout twice more | typed errors, exponential backoff with jitter, `Retry-After`, no retry on permanent failures |
| No response validation | an HTML error page delivered with HTTP 200 was written out as a `.pdb` file | content validated before acceptance |
| No caching | repeating an analysis re-hit the public services | checksum-verified on-disk cache |
| Biopython required | installation failures on machines without a build toolchain | no required dependencies; Biopython optional, used for cross-validation |

---

## 4. Scientific findings

Reviewing the calculations against their published sources found seven
defects. Each produced plausible output, which is why none had been noticed.
They are listed with their fixes in [`../CHANGELOG.md`](../CHANGELOG.md); the
two most consequential:

**Reading-frame shift.** Ambiguity codes were stripped from the sequence
before translation. Removing a single `N` shifts every subsequent codon, so
`ATGNAAAGGTTT` produced a peptide that looked like a normal translation and
was wrong from residue two onward. Because the output was a plausible protein
string, nothing flagged it. Ambiguity codes now stay in place and their
codons translate to `X`.

**Melting temperature.** The Wallace rule was applied to sequences of any
length without comment. It is valid for oligonucleotides of roughly 14–20 nt;
applied to a gene it returns a large number with no physical meaning. v2.0
uses SantaLucia (1998) nearest-neighbour thermodynamics with salt correction
in the appropriate range, reports which method and conditions produced the
number, and attaches a warning outside it.

Two further points concern presentation rather than arithmetic:

* The secondary-structure display was described in wording that invited it to
  be read as a prediction. It is a propensity sketch. It is now labelled as
  such in the interface, in the exported data and in the documentation.
* pLDDT and crystallographic B-factors were distinguished by the numeric
  range of the B-factor column, which cannot tell a well-ordered crystal
  structure from a predicted model. Provenance now decides the label.

---

## 5. Evidence findings

| Claim in the v1.0 manuscript | What supported it | What v2.0 provides |
|---|---|---|
| "25/25 structures retrieved" | a count of HTTP 200 responses | mapping-level distribution with identity and coverage per protein; the 7 accession-discordant cases are adjudicated individually |
| "robust error handling" | the presence of `try`/`except` | 13-scenario fault-injection experiment over 1 300 requests, reproducible offline in seconds |
| "9/25 ESM Atlas requests returned 504" | an observation | the same failure modes injected deterministically, with a measured 80 % recovery rate |
| "integrates multiple analyses" | a feature list | interaction counts against a published manual protocol |
| "user-friendly" | assertion | a pre-registered usability protocol; **no results claimed until it is run** |

The last row is the one to be careful about. v2.0 does not claim a usability
result. It publishes the protocol and confines the present claims to
interaction counts and machine timings, which need no participants.

---

## 6. What v1.0 got right

Worth recording, since the rewrite kept it:

* The problem is real. Fragmented sequence-analysis workflows with an unchecked
  structure lookup at the end are how a great deal of routine work is done.
* Querying three structural resources with fallback is the right shape for a
  client, because their coverage genuinely differs.
* A desktop GUI is the right delivery mechanism for users who will not write
  a Biopython script.
* The four-tab organisation was sound enough to keep.

---

## 7. How to verify this audit

The v1.0 source is at the `v1.0-archive` tag. The mis-mapping claims can be
checked directly:

```bash
# What v2.0 reports for the preproinsulin case
bioseqinsight structure P01308 --all --json
```

The returned accessions and measured identity are in the output. Any protein
in `benchmarks/data/benchmark_candidates.csv` can be checked the same way,
and the full live benchmark regenerates the table:

```bash
python benchmarks/scripts/run_structure_benchmark.py --out results/audit_check
```

Dispute is welcome and should take the form of a run whose output disagrees
with what is written here.
