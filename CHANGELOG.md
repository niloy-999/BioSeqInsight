# Changelog

All notable changes to BioSeqInsight are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project
follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.0.2] — 2026-09-23

Two more issues found from a second independent live run of the structure
benchmark, on top of the four fixed in 2.0.1 the day before.

### Fixed

1. **A protein with no structure from any provider vanished from the
   mapping-level distribution instead of being counted as M0.**
   `run_structure_benchmark.py`'s `best_mapping_distribution` was built by
   starting each protein at an implicit "M0" and only ever recording an
   entry when some provider beat it — so a protein that never beat M0
   anywhere was never written into the results dictionary at all, and a
   value that is never written is not counted by anything that only
   iterates the dictionary's values. A second live run over the same
   69-protein set produced a table reading `M0:0, M1:0, M2:0, M3:1, M4:67`
   — summing to 68 against 69 resolved candidates, with the one protein
   that genuinely got nothing simply absent rather than tallied under M0.
   The two headline rates (`exact_mapping_rate_percent`,
   `any_structure_rate_percent`) were unaffected, because neither is derived
   from `distribution["M0"]` — but the table itself, which is exactly the
   kind of thing meant to go in a manuscript, did not reconcile. Every
   resolved accession is now seeded into the results dictionary as M0
   before any row is processed, and the script asserts the distribution
   sums to the resolved count before returning, so this class of defect
   cannot recur silently again. (`benchmarks/scripts/run_structure_benchmark.py`)
2. **`success_rate_percent` conflated three different things under one
   number.** For a provider like ESM Atlas, which refuses sequences over a
   length limit before making any request, `success_rate_percent` blended
   genuine service failures, genuine coverage gaps, and sequences that were
   never attempted at all into a single denominator — a live run reported
   "66.7%" for ESM Atlas in a way that reads as unreliability but was
   actually 100% success among the 46 sequences short enough to try, diluted
   by 23 the software correctly never attempted. Added
   `attempted_requests` and `attempted_success_rate_percent`, computed over
   only the sequences actually sent to the service, reported alongside the
   original (unchanged) `success_rate_percent` rather than replacing it —
   the two answer different questions and a report should be explicit about
   which one it means. (`benchmarks/scripts/run_structure_benchmark.py`)

Both were found by checking that a summary table's own numbers reconciled
with each other, not by anything the offline test suite could exercise —
the offline suite's synthetic scenarios happened never to include a protein
that scored M0 against every provider, so the dictionary-seeding gap was
never triggered by any existing test. Six regression tests
(`tests/test_benchmark_scripts.py::TestBestMappingDistribution`,
`::TestAttemptedSuccessRate`) now cover both directly, including the exact
"one candidate genuinely gets nothing" shape that a live run hit.

---

## [2.0.1] — 2026-09-23

Four issues found by running the software against real network services and
a real Biopython installation, none of which the offline test suite alone
could catch. All four are fixed and confirmed by a second live run.

### Fixed

1. **AlphaFold DB retrieval failed for every accession.** The provider
   constructed the model file URL by guessing a version number
   (`AF-{accession}-F1-model_v4.pdb`, falling back to `v3`). Once EBI
   republished AlphaFold DB past those versions, every request 404'd — the
   first live run against 69 real proteins showed **0% AlphaFold
   retrieval**, with 60 of 69 reported as "not found," indistinguishable
   from a genuine coverage gap unless something checks the authoritative
   source. The provider now queries the public prediction API
   (`/api/prediction/{accession}`) first and fetches whichever model file URL
   it returns — AlphaFold DB turned out to already be on `v6` by the time of
   the second live run, confirming the version-guessing approach could never
   have self-corrected. (`structures/alphafold.py`)
2. **A clean two-call success was misreported as a recovery from failure.**
   `StructureResult.recovered` was defined as `ok and len(attempts) > 1`,
   safe only because every existing provider made exactly one logical HTTP
   request. The AlphaFold fix above introduced a provider that legitimately
   makes two calls (discover the model URL, then download it) even on a
   completely clean run, and the second live run's first pass showed **67 of
   67 successful AlphaFold retrievals flagged "recovered after failure,"**
   none of which had actually failed. `recovered` is now computed from
   whether any individual attempt actually errored or returned a 4xx/5xx
   status. (`models/results.py::StructureResult.recovered`)
3. **Isoelectric point test overclaimed agreement with Biopython.** The
   parity test asserted agreement within 0.05 pH units with a comment
   claiming the two implementations "implement the same charge model."
   Running against a real Biopython installation showed differences of
   0.23-0.26 pH units on two of four test sequences: the pKa tables genuinely
   differ, which is expected and documented in the pI-calculator literature,
   but the test's claim was wrong. The test now checks agreement within a
   defensible tolerance (0.5 pH units) and adds a self-consistency check —
   net charge is zero at the reported pI — which does not depend on
   agreeing with any particular external tool's table. (`core/protein.py`,
   `tests/test_biopython_parity.py`)
4. **The structure-benchmark summary silently dropped a whole status
   category.** `run_structure_benchmark.py`'s per-provider summary reported
   `retrieved`, `not_found`, `timeouts` and `service_errors`, which do not
   sum to `requests` when a provider refuses input before making any HTTP
   call. ESM Atlas does exactly this for sequences over its length limit
   (`esm_max_length`, 400 aa by default): the second live run showed 46 of 69
   ESM Atlas requests "retrieved" with the other 23 unaccounted for anywhere
   in the output. The 23 are exactly the candidates longer than 400 aa. The
   summary now reports `invalid_input` and `skipped` as their own fields and
   an `accounting_check` that states outright whether every request in the
   run landed in a counted category, so this cannot recur silently.
   (`benchmarks/scripts/run_structure_benchmark.py`)

**Confirmed by a second live run** (`benchmarks/results/structure_benchmark_2026-09-23_summary.json`)
after all four fixes: AlphaFold retrieval **0% → 98.6%** (68/69), overall
exact (M4) mapping rate **98.6%**, any-structure rate **100%** across 69
proteins, `recovered_after_failure` a plausible `0` across all three
providers, and every provider's accounting reconciling exactly. Both earlier
runs are kept in `benchmarks/results/` as an audit trail; see
`benchmarks/results/README.md` for the full before/after sequence.

None of these four were reachable from the offline test suite on their own:
each depended on either real infrastructure having moved past a hard-coded
assumption, or on a real third-party library being installed, or on running
the full candidate set rather than a handful of unit-test fixtures. This is
the argument for `docs/benchmarking.md`'s insistence on a dated live run
before any retrieval figure goes into the manuscript, made concrete.

---

## [2.0.0] — 2026-09-21

A rewrite. The scientific calculations were re-derived and re-tested rather
than moved, the architecture was rebuilt in layers, and the central new
capability — validation of sequence-to-structure mapping — did not exist in
any form in version 1.0.

Version 1.0 is preserved at the `v1.0-archive` tag so that results published
from it can still be reproduced.

### Added

**Structure validation (the reason for this release)**

- **M0–M4 sequence-to-structure mapping taxonomy.** Every retrieved structure
  is aligned against the query sequence and cross-checked against its
  database accession. The outcome is reported as one of five levels in the
  GUI, the CLI and every exported table. Only M4 (accession confirmed *and*
  identity ≥ 99 % at ≥ 95 % coverage) supports the claim that a structure is
  the queried protein.
- Needleman–Wunsch and semi-global pairwise alignment with affine gap
  penalties (`core/alignment.py`), used for identity and coverage
  measurement. Semi-global alignment means a PDB chain that is a perfect
  fragment of a longer UniProt sequence scores 100 % identity with reduced
  coverage, which is the biologically correct reading.
- Chain-level parsing of coordinate files with SEQRES/ATOM reconciliation and
  modified-residue mapping (`structures/pdbio.py`).
- Explicit warnings when a retrieved structure is a homologue, an orthologue,
  a fragment or a different molecule.

**Architecture**

- Layered package under `src/`: `core`, `models`, `services`, `structures`,
  `workflows`, `io`, `gui`, `cli`. Dependencies point one way.
- `core` has no I/O, no network and no GUI imports, and uses only the
  standard library.
- Typed, JSON-serialisable result dataclasses (`models/results.py`) replace
  the preformatted display strings that v1.0 returned from every function.
- Pluggable HTTP transports (`requests`, `urllib`, and a scripted fake for
  tests), selected at run time.

**Interfaces**

- A complete command-line interface with eight subcommands and meaningful
  exit codes. Everything the GUI can do is reachable without a display.
- Batch analysis: parallel, fault-isolated, with CSV, TSV, JSON and
  self-contained HTML export.
- Projects: a directory holding inputs, retrieved structures, results, logs
  and a snapshot of the software version and settings, exportable as a zip
  that re-imports into a working project.
- GUI gained a Batch tab, a mapping-level badge, per-resource attempt
  reporting, and background execution so the window no longer freezes during
  network calls.

**Reliability**

- Typed error hierarchy replacing bare exceptions.
- Exponential backoff with jitter, a configurable attempt budget, and
  `Retry-After` support. Permanent failures (404, 403) are not retried.
- Content-level response validation: a gateway answering HTTP 200 with an
  HTML error page is treated as a failure instead of being written out as a
  structure file.
- Multi-resource fallback: RCSB → AlphaFold DB → ESM Atlas, with early exit
  on an exact match only.
- Checksum-verified on-disk response cache with TTL and CLI management.
- Structured logging (optionally JSON) that records a salted hash and length
  of a sequence rather than the sequence itself.

**Configuration**

- `Settings` with a documented precedence chain: defaults → file →
  environment → command line, validated on load and recorded in every project
  and JSON export.

**Testing and CI**

- 432 tests, 89.6 % statement coverage (90.1 % across `core`, `structures`
  and `services`). No test touches the network.
- `tests/test_biopython_parity.py` cross-validates every dependency-free
  calculation against Biopython; it skips when Biopython is absent and runs
  in a dedicated CI job.
- GitHub Actions: 12-way OS × Python matrix, a no-dependency job that fails
  if a third-party import creeps into the core, a Biopython cross-validation
  job, a requests-transport job, an offline benchmark job, and a job that
  installs the built wheel into a clean environment outside the source tree.

**Benchmarks**

- Deterministic offline fault-injection harness covering 13 failure
  scenarios.
- Offline performance and scalability benchmark with a fixed seed.
- Workflow comparison deriving interaction counts from a published manual
  protocol.
- Live structure benchmark over a curated candidate set, which resolves
  lengths from UniProt at run time rather than trusting a hand-typed column.
- Pre-registered usability study protocol (protocol only; no participant data
  is included, and none is claimed).

**Documentation**

- README, installation guide, user guide, six tutorials, architecture note,
  developer guide, benchmarking guide, mapping-level reference, contribution
  guide, security policy, and an audit of version 1.0.

### Changed

- **Biopython is no longer required.** Every sequence and protein calculation
  is implemented against the standard library, so the package installs and
  its full test suite runs anywhere Python runs. Biopython becomes an
  optional extra used to cross-validate those implementations — a stronger
  position than depending on it, because agreement is now demonstrated rather
  than assumed.
- Melting temperature now uses SantaLucia (1998) nearest-neighbour
  thermodynamics with salt correction for 14–60 nt duplexes, and reports the
  method and conditions used. v1.0 applied the Wallace rule to sequences of
  any length without comment.
- The secondary-structure display is now labelled a *propensity sketch*
  everywhere it appears, including in the exported data. v1.0's wording
  invited it to be read as a prediction.
- pLDDT and crystallographic B-factors are distinguished by provider
  provenance rather than by guessing from the numeric range, and are labelled
  differently in the interface.
- Retrieval reports every resource queried, with attempt counts, rather than
  a single success or failure line.
- Imports are package-relative throughout; v1.0 could only run from its own
  directory.

### Fixed

Defects found by writing the test suite. Each was present in v1.0 and each
produced plausible-looking wrong output rather than an error, which is why
none had been noticed.

1. **Reading-frame shift from ambiguity-code stripping.** Removing `N` before
   translation shifted every downstream codon. `ATGNAAAGGTTT` translated to
   `MKRF`-like nonsense instead of `MXRF`. Ambiguity codes now stay in place
   and their codons translate to `X`. This affected translation, ORF
   detection and motif coordinates.
2. **Silent deletion of invalid characters.** A stray symbol such as `@` was
   stripped and the sequence analysed anyway. Genuine symbols are now
   reported; digits, gaps and whitespace, which legitimately appear in text
   copied from database records, are still ignored.
3. **Wrong alphabet for very short input.** Two- and three-letter strings
   were confidently classified as protein. Since every amino-acid letter is
   also a possible nucleotide ambiguity code, inputs under four letters are
   now reported as undetermined instead of guessed.
4. **Accession parsing of FASTA headers.** `sp|P24941|CDK2_HUMAN` yielded
   `CDK2_HUMAN`, so accession cross-checks silently failed. The parser now
   finds the accession-shaped field wherever it sits.
5. **pLDDT and B-factor conflation.** Confidence was inferred from the
   numeric range of the B-factor column, which cannot distinguish a
   well-ordered crystal structure from a predicted model. Provenance now
   decides the label.
6. **Measured 0 % identity promoted to a "related record".** A structure that
   had been compared and found unrelated was classified M2. It is now M1 with
   an explicit warning; M2 means related, and must not be used for something
   that is not.
7. **Viewer breakage on coordinate files containing backticks or
   backslashes.** The v1.0 viewer interpolated PDB text into a JavaScript
   template literal. Coordinates are now embedded with `json.dumps`.

### Security

- Zip extraction rejects entries that would escape the destination directory,
  with a regression test.
- Generated HTML escapes user-controlled identifiers.
- The stdlib transport refuses non-HTTPS URLs.
- Cached payloads are verified against a stored SHA-256 digest before use.

### Known limitations

Documented rather than worked around; see README section 14 for the full
list. The software does not predict structure, the propensity sketch is not a
prediction, alignment is for identity verification rather than homology
search, and the human-factors claims rest on interaction counts and machine
timings, not on a completed usability study.

---

## [1.0.0] — 2026

Initial release: a Tkinter application in three modules providing GC content,
reverse complement, transcription, translation, ORF detection, Wallace-rule
melting temperature, protein physicochemical properties, a secondary-structure
sketch, and structure retrieval from RCSB PDB, AlphaFold DB and ESM Atlas.

Retained for reproducibility at the `v1.0-archive` tag. Its retrieval results
should be read alongside [`docs/SOFTWARE_AUDIT.md`](docs/SOFTWARE_AUDIT.md),
which shows that 7 of its 25 reported successful retrievals were a different
molecule from the one queried.
