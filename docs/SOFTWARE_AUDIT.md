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

**The v1.0 manuscript reported 25 of 25 proteins successfully mapped to
RCSB PDB entries. Re-examining those 25 mappings shows that 7 of them
returned a different molecule from the one queried.**

The software had no way to notice. It performed a sequence search, took a
hit, downloaded the coordinates, and reported success. No step compared the
deposited chain with the query, and no step checked whether the entry's
UniProt cross-reference matched the accession that had been asked about.

The clearest case: querying human **preproinsulin (P01308)** returned a
structure of **insulin-degrading enzyme (P14735)** — the enzyme that destroys
the hormone rather than the hormone itself. Both are legitimate sequence-search
neighbours. Only one is the answer to the question.

Other categories of mis-mapping in the same set:

| Category | Example | Why the search returned it |
|---|---|---|
| Substrate returned instead of the protein | preproinsulin → insulin-degrading enzyme | local sequence similarity in the bound-peptide region |
| Complex partner instead of the subunit | GroES (P0A6F9) → GroEL | co-deposited in the same entry |
| Non-human orthologue for a human query | β2-microglobulin (P61769) → non-human entry | high cross-species conservation |
| Fragment presented as the full protein | multi-domain queries → single-domain entries | the deposited construct is a fragment |

None of these are failures of the RCSB search, which did what it was asked.
They are failures of a client that treated "a structure came back" as "the
right structure came back".

**What v2.0 does about it.** Every retrieval is aligned against the query and
cross-checked against the accession, and the result is reported on the M0–M4
scale. The preproinsulin case now comes back as M1 with an explicit warning
naming the accession actually returned, instead of as a success.

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
| "25/25 structures retrieved" | a count of HTTP 200 responses | mapping-level distribution with identity and coverage per protein; the 7 mis-mappings are visible |
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
