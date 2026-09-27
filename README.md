# BioSeqInsight 2.0

**An integrated desktop workflow for DNA sequence analysis and *validated* protein structure-resource retrieval.**

[![tests](https://github.com/niloy-999/BioSeqInsight/actions/workflows/tests.yml/badge.svg)](https://github.com/niloy-999/BioSeqInsight/actions/workflows/tests.yml)
[![lint](https://github.com/niloy-999/BioSeqInsight/actions/workflows/lint.yml/badge.svg)](https://github.com/niloy-999/BioSeqInsight/actions/workflows/lint.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

BioSeqInsight combines routine nucleotide and protein analysis with retrieval
from RCSB PDB, the AlphaFold Protein Structure Database and the ESM Atlas
fold service — and, unlike a browser tab, it tells you whether the structure
it retrieved is actually the protein you asked for.

```bash
pip install bioseqinsight
bioseqinsight structure P0CG48
```

```
[RCSB] ok
  identifier     1UBQ
  mapping        M4 Exact accession and sequence match
  identity       100.0%   coverage 100.0%   chain A

BEST RESULT      rcsb 1UBQ
MATCH STATUS     M4 - Exact accession and sequence match
Exact identity    YES
```

---

## Contents

1. [Why this exists](#1-why-this-exists)
2. [What is new in 2.0](#2-what-is-new-in-20)
3. [Installation](#3-installation)
4. [Quick start](#4-quick-start)
5. [The mapping levels (M0–M4)](#5-the-mapping-levels-m0m4)
6. [Command-line interface](#6-command-line-interface)
7. [Graphical interface](#7-graphical-interface)
8. [Batch analysis and projects](#8-batch-analysis-and-projects)
9. [Configuration](#9-configuration)
10. [Architecture](#10-architecture)
11. [Testing](#11-testing)
12. [Benchmarks and reproducibility](#12-benchmarks-and-reproducibility)
13. [Fault tolerance](#13-fault-tolerance)
14. [Scope and limitations](#14-scope-and-limitations)
15. [Comparison with existing tools](#15-comparison-with-existing-tools)
16. [Documentation](#16-documentation)
17. [Contributing](#17-contributing)
18. [Citation](#18-citation)
19. [Licence and acknowledgements](#19-licence-and-acknowledgements)

---

## 1. Why this exists

Characterising a gene and finding a structure for its product is routine
work, and it is routinely done by hand across four or five websites: a GC
calculator, a Tm calculator, an ORF finder, a protein-parameter page, then
UniProt and the PDB. The sequence gets pasted five times. The parameters each
tool defaulted to are not recorded anywhere. And the step that establishes
whether the retrieved structure is your protein — aligning the deposited
chain against your query — is the easiest one to skip.

Skipping it is not hypothetical. The audit in [`docs/SOFTWARE_AUDIT.md`](docs/SOFTWARE_AUDIT.md)
re-examined the 25-protein benchmark published for version 1.0 of this
software, which reported 25/25 successful structure retrievals. Seven of
those 25 were a different molecule. Querying human preproinsulin returned
insulin-*degrading* enzyme, and the software reported success.

BioSeqInsight 2.0 is built around that finding. Every retrieved structure is
aligned against the query sequence and cross-checked against the accession,
and the result is reported on an explicit five-level scale. The software will
tell you it failed rather than hand you the wrong protein.

---

## 2. What is new in 2.0

Version 2.0 is a rewrite, not an increment.

| | v1.0 | v2.0 |
|---|---|---|
| **Structure validation** | none | M0–M4 taxonomy: alignment identity, coverage and accession cross-check |
| **Architecture** | 3 flat modules, 1 242 lines | layered package, ~7 000 lines, GUI-free core |
| **Return values** | preformatted display strings | typed dataclasses, JSON-serialisable |
| **Tests** | 108 lines, no framework | 432 tests, 89.6 % statement coverage |
| **CI** | none | 12-way OS × Python matrix, plus lint, packaging and benchmark jobs |
| **Interfaces** | Tkinter only | Tkinter **and** a complete CLI |
| **Batch processing** | none | parallel, fault-isolated, CSV/TSV/JSON/HTML export |
| **Reproducibility** | none | projects store inputs, settings and software version; re-importable |
| **Fault tolerance** | widen the timeout, try 3× | typed errors, exponential backoff with jitter, `Retry-After`, payload validation, multi-resource fallback |
| **Dependencies** | Biopython required | **none required**; Biopython optional, used for cross-validation |
| **Caching** | none | checksum-verified on-disk cache with TTL |

Version 1.0 remains available on the `v1.0-archive` tag so the published
results can still be reproduced.

---

## 3. Installation

Python 3.10 or newer. No compiler, no scientific stack.

```bash
pip install bioseqinsight
```

From source:

```bash
git clone https://github.com/niloy-999/BioSeqInsight.git
cd BioSeqInsight
pip install -e .
```

Optional extras:

```bash
pip install "bioseqinsight[http]"        # use requests instead of urllib
pip install "bioseqinsight[validation]"  # Biopython, for the parity tests
pip install "bioseqinsight[dev]"         # pytest, ruff, black, mypy
```

The GUI needs Tkinter, which ships with most Python builds. On Debian and
Ubuntu: `sudo apt install python3-tk`. Everything the GUI does is also
available from the command line, so a headless server is fully supported.

Verify:

```bash
bioseqinsight version
python -m pytest -q
```

Full notes, including Windows and macOS specifics: [`docs/installation.md`](docs/installation.md).

---

## 4. Quick start

```bash
# Nucleotide analysis: GC, Tm, ORFs, composition
bioseqinsight dna --file examples/dna/demo_gene.fasta

# Protein descriptors
bioseqinsight protein --file examples/proteins/ubiquitin.fasta --sketch

# Retrieve a structure and validate its identity
bioseqinsight structure P0CG48

# A whole FASTA file to a table
bioseqinsight batch examples/proteins/benchmark_subset.fasta \
    --structures --formats csv,html

# Launch the desktop interface
bioseqinsight-gui
```

As a library:

```python
from bioseqinsight import analyze_sequence, ProteinAnalyzer, StructureManager

result = analyze_sequence("ATGGCTAGCAAAGGTTTCCCG...", identifier="demo")
print(result.gc_percent, len(result.orfs))

outcome = StructureManager().retrieve("P24941")
best = outcome.best
print(best.identifier, outcome.best_mapping.value, best.validation.identity_percent)
```

Every result object has `.to_dict()` and serialises to JSON, so results can
be piped onward rather than screen-scraped.

---

## 5. The mapping levels (M0–M4)

The central claim of this software is that it never reports a structure as
your protein without evidence. That evidence is summarised in one level,
shown in the GUI, the CLI and every exported table.

| Level | Meaning | Safe to state "this is the structure of my protein"? |
|---|---|---|
| **M4** | Exact: the accession is cross-referenced by the entry **and** identity ≥ 99 % with coverage ≥ 95 % | **Yes** |
| **M3** | High confidence: identity ≥ 95 %, coverage ≥ 90 %, accession not confirmed | Probably, state the evidence |
| **M2** | Related record: a homologue, orthologue or fragment | **No** |
| **M1** | Retrieved but identity not established, or measured and clearly different | **No** |
| **M0** | No structural record returned | — |

Only M4 justifies the unqualified claim. The thresholds are configurable and
recorded with every result, so a reader can see which ones produced a number.

Full rationale, worked examples and the threshold sensitivity analysis:
[`docs/mapping-levels.md`](docs/mapping-levels.md).

---

## 6. Command-line interface

The CLI is not a convenience wrapper; it is how the published evaluation is
reproduced, on a machine with no display.

| Command | Purpose |
|---|---|
| `dna` | nucleotide analysis: GC, Tm, ORFs, motifs, composition |
| `protein` | molecular weight, pI, GRAVY, aromaticity, extinction coefficient |
| `structure` | retrieve and validate from RCSB, AlphaFold DB, ESM Atlas |
| `batch` | run a FASTA file, export CSV/TSV/JSON/HTML |
| `project` | create, inspect, export and import reproducible projects |
| `cache` | inspect or clear the response cache |
| `config` | show or save the effective settings |
| `version` | version, platform and optional-dependency report |

Exit codes: `0` success, `2` user error, `3` runtime failure — so it composes
in shell scripts and CI.

```bash
bioseqinsight dna ATGGCTAGC... --motif GAATTC --min-orf 50 --json
bioseqinsight structure P24941 --providers rcsb,alphafold --all --json
bioseqinsight batch seqs.fasta --structures --workers 8 --project ./run01
```

---

## 7. Graphical interface

Four tabs — Sequence Analysis, Protein Analysis, Structure Retrieval, Batch
Analysis — over the same core API the CLI uses. Network calls run on a
background thread, so the window never freezes while a service is slow.

The Structure tab shows every resource that was queried, how many attempts
each took, whether it recovered from a failure, and the mapping level with
the identity and coverage behind it. Retrieved coordinates open in an
embedded 3Dmol.js viewer or an external program such as PyMOL or ChimeraX.

Walkthrough with screenshots: [`docs/user-guide.md`](docs/user-guide.md).

---

## 8. Batch analysis and projects

Batch mode processes a FASTA file in parallel. One malformed record does not
stop the run; it is reported in the output table with its error.

A **project** is a directory that makes a run repeatable:

```
project/
├── project.json      manifest: name, record count, notes, format version
├── metadata.json     software version + complete settings snapshot
├── sequences/        the exact inputs analysed
├── structures/       retrieved coordinates
├── results/          exported tables
└── logs/             structured run logs
```

`bioseqinsight project export` produces a zip that re-imports into a working
project, which makes "here is everything needed to check my analysis" a
single action rather than a folder-assembly exercise.

---

## 9. Configuration

Settings resolve in order: **defaults → file → environment → command line**,
and the effective values are recorded in every project and every JSON export.

```bash
bioseqinsight config                                    # show effective settings
bioseqinsight config --save ~/.bioseqinsight/settings.json
export BIOSEQINSIGHT_MAX_ATTEMPTS=5
bioseqinsight --offline --no-cache structure P24941
```

Configurable: endpoints, timeouts, retry budget and backoff, cache location
and TTL, download directory, provider order, mapping thresholds, Tm
conditions, genetic-code table, batch workers, logging.

---

## 10. Architecture

```
bioseqinsight/
├── core/         sequence and protein computation — no I/O, no network, no GUI
├── models/       typed result dataclasses shared by every layer
├── services/     HTTP transport, retry policy, cache, logging
├── structures/   providers (RCSB, AlphaFold, ESM Atlas, UniProt) + validation + manager
├── workflows/    batch runner, project management
├── io/           CSV, TSV, JSON and self-contained HTML export
├── gui/          Tkinter views (formatting logic kept Tkinter-free, so it is testable)
└── cli.py        command-line entry point
```

Dependencies point one way: `gui` and `cli` → `workflows` → `structures` →
`services` → `core`. `core` imports nothing from the layers above it and
nothing from outside the standard library, which is why the scientific code
can be tested exhaustively offline.

See [`docs/architecture.md`](docs/architecture.md) for the reasoning behind
each boundary.

---

## 11. Testing

```bash
python -m pytest                              # 432 tests, offline, ~20 s
python -m pytest --cov=bioseqinsight           # with coverage
python -m pytest tests/test_validation.py -v   # the mapping taxonomy
```

**432 tests, 89.6 % statement coverage** (90.1 % across `core`, `structures`
and `services`). No test contacts the network: HTTP is exercised through a
scripted transport, which makes failure-path testing deterministic instead of
dependent on a public service failing at the right moment.

The `tests/test_biopython_parity.py` module cross-validates every
dependency-free calculation against Biopython. It skips when Biopython is
absent and runs in a dedicated CI job, so the comparison always happens
before a release.

Writing the suite found seven genuine defects in code that appeared to work,
including a reading-frame shift caused by stripping ambiguity codes before
translation. They are listed in [`CHANGELOG.md`](CHANGELOG.md).

---

## 12. Benchmarks and reproducibility

Everything in `benchmarks/` regenerates from a single command.

| Script | Network? | Measures |
|---|---|---|
| `run_fault_injection.py` | no | recovery from 13 scripted failure modes |
| `run_performance_benchmark.py` | no | throughput, scaling, alignment cost, cache latency |
| `run_workflow_comparison.py` | no | interaction cost versus the documented manual workflow |
| `run_structure_benchmark.py` | **yes** | live retrieval and mapping levels across the curated protein set |

```bash
python benchmarks/scripts/run_fault_injection.py
python benchmarks/scripts/run_performance_benchmark.py
python benchmarks/scripts/run_structure_benchmark.py --out results/run_$(date +%F)
```

The candidate protein set (`benchmarks/data/benchmark_candidates.csv`) stores
accessions and rationale but **no sequence lengths**: lengths are fetched from
UniProt at run time and the size strata are computed from them. A published
table can therefore never silently disagree with the sequences it describes.

Measured offline on the reference machine (Python 3.12, x86-64):

* **850–900 sequences/second** for local analysis; per-sequence cost varies by
  a factor of 1.18 from 10 to 1 000 sequences, i.e. throughput is linear.
* **0.2 ms** cache read for a 200 kB coordinate file.
* **25 ms** from process start to first completed analysis.

Details and how to re-run: [`docs/benchmarking.md`](docs/benchmarking.md).

---

## 13. Fault tolerance

Version 1.0 reported that 9 of 25 ESM Atlas requests returned HTTP 504.
Reporting an outage is not an engineering result; handling it is.

Failures are typed, retried with exponential backoff and jitter, and bounded
by a configurable budget. `Retry-After` is honoured. Permanent failures (404,
403) are not retried. Responses are validated as *content*, so a gateway that
answers HTTP 200 with an HTML error page is treated as the failure it is. If
one resource cannot help, the manager falls back to the next.

The deterministic fault-injection experiment covers 13 scenarios across 1 300
requests. Every scenario behaved as specified, and **800 of the 1 000 requests
that began with a failure recovered (80.0 %)** — the remaining 200 being the
scenarios designed to be unrecoverable.

| Scenario | Outcome |
|---|---|
| Single 504, then success | 100/100 recovered |
| Two 504s, then success | 100/100 recovered |
| 500 / 502 / 503, then success | 100/100 recovered each |
| Rate limited (429 + `Retry-After`) | 100/100 recovered |
| Read timeout, then success | 100/100 recovered |
| HTTP 200 carrying an HTML error page | 100/100 detected and recovered |
| Persistent 504 | 100/100 correctly unrecoverable at 3 attempts |
| 404 / 403 | not retried (1 attempt), as intended |

Reproduce in about two seconds, with no network:

```bash
python benchmarks/scripts/run_fault_injection.py
```

---

## 14. Scope and limitations

Stated plainly, because a tool that overstates itself is worse than one that
does less.

**BioSeqInsight does not predict protein structure.** It retrieves models
from AlphaFold DB and submits sequences to the ESM Atlas fold service. The
prediction is theirs; cite them, not this software.

**The secondary-structure display is a propensity sketch, not a prediction.**
It applies Chou–Fasman-style propensities for illustration. It is labelled as
such in the interface and must not be reported as a prediction. Use DSSP on a
real structure, or a dedicated predictor.

**Other limits.**

* Sequence search reflects what RCSB's API returns; BioSeqInsight validates
  the result but does not re-rank the search.
* ESM Atlas enforces a sequence-length limit; longer sequences are refused
  before the request rather than failing remotely.
* Melting temperature uses SantaLucia (1998) nearest-neighbour parameters for
  duplexes of 14–60 nt. Outside that range a warning is attached to the
  result, and long sequences do not obey a two-state model.
* Alignment is Needleman–Wunsch / semi-global with an affine gap penalty,
  capped at 4 million cells. It is for identity verification, not homology
  searching; use BLAST or HMMER for that.
* Isoelectric point uses a fixed pKa set and ignores local structural
  effects, as all sequence-based pI calculators do.
* Network operations depend on third-party services whose availability and
  APIs can change.
* The human-factors claims rest on interaction counts and machine timings
  only. The usability study
  ([`benchmarks/protocols/usability_study_protocol.md`](benchmarks/protocols/usability_study_protocol.md))
  is specified but its results are not included here, and nothing in this
  repository claims otherwise.

---

## 15. Comparison with existing tools

| | Biopython | EMBOSS | Web tools | UCSF ChimeraX | **BioSeqInsight** |
|---|---|---|---|---|---|
| Sequence analysis | ✔ library | ✔ CLI | ✔ per tool | partial | ✔ GUI + CLI + library |
| Structure retrieval | partial | ✖ | manual | ✔ | ✔ multi-resource with fallback |
| **Identity validation of the retrieved structure** | do it yourself | do it yourself | ✖ | ✖ | **✔ automatic, M0–M4** |
| Batch processing | write a script | ✔ | ✖ | scriptable | ✔ built in |
| Records the settings used | ✖ | ✖ | ✖ | session file | ✔ project snapshot |
| Usable without programming | ✖ | partial | ✔ | ✔ | ✔ |
| Installation weight | library | system package | none | large | `pip install`, no dependencies |

BioSeqInsight is not a replacement for Biopython, and a bioinformatician
scripting a pipeline should use Biopython. The gap it fills is the user who
will not write that script and who currently does this work by hand across
browser tabs, with no record of what was done and no check on whether the
structure is the right molecule.

---

## 16. Documentation

| Document | Contents |
|---|---|
| [`docs/installation.md`](docs/installation.md) | platform-specific installation and troubleshooting |
| [`docs/user-guide.md`](docs/user-guide.md) | GUI and CLI walkthrough |
| [`docs/tutorials/`](docs/tutorials/) | six worked tutorials |
| [`docs/mapping-levels.md`](docs/mapping-levels.md) | the M0–M4 taxonomy in full |
| [`docs/architecture.md`](docs/architecture.md) | layer boundaries and why they are where they are |
| [`docs/developer-guide.md`](docs/developer-guide.md) | extending the software, adding a provider |
| [`docs/benchmarking.md`](docs/benchmarking.md) | reproducing every published number |
| [`docs/SOFTWARE_AUDIT.md`](docs/SOFTWARE_AUDIT.md) | the v1.0 audit, including the 7/25 mis-mappings |
| [`CHANGELOG.md`](CHANGELOG.md) | every change, including the defects the tests found |

---

## 17. Contributing

Contributions are welcome. Please read
[`CONTRIBUTING.md`](CONTRIBUTING.md); in short: new behaviour needs a test,
the suite must stay offline and deterministic, and public functions need
docstrings that say what the function does *and* what it does not.

Security issues: [`SECURITY.md`](SECURITY.md).

---

## 18. Citation

If BioSeqInsight is useful in your work, please cite it — and please cite the
structural resources separately, because the science in a retrieved structure
is theirs:

* **RCSB PDB** — Berman et al. (2000) *Nucleic Acids Research* 28:235–242
* **AlphaFold DB** — Varadi et al. (2022) *Nucleic Acids Research* 50:D439–D444
* **AlphaFold** — Jumper et al. (2021) *Nature* 596:583–589
* **ESMFold / ESM Atlas** — Lin et al. (2023) *Science* 379:1123–1130
* **UniProt** — The UniProt Consortium (2023) *Nucleic Acids Research* 51:D523–D531

Machine-readable metadata for this software is in
[`CITATION.cff`](CITATION.cff).

---

## 19. Licence and acknowledgements

MIT — see [`LICENSE`](LICENSE).

BioSeqInsight is a client for public scientific infrastructure maintained by
the RCSB PDB, EMBL-EBI and DeepMind (AlphaFold DB), Meta AI (ESM Atlas) and
the UniProt Consortium. It implements none of their science and redistributes
none of their data. Please observe each resource's terms of use and keep
automated request rates reasonable.

Structure viewing uses [3Dmol.js](https://3dmol.csb.pitt.edu/) (Rego &
Koes, 2015), loaded from a CDN in the generated viewer page.
