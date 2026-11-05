# Validation benchmark protocol v1

**Status: frozen before Experiment A is run.** Changes after freezing are
logged in the Amendments section at the bottom, each with a reason. Do not
edit the criteria above that line to make a disappointing result look
better; if a genuine methodological problem is found, document it here,
re-run, and report both runs.

This protocol implements the six-experiment plan: a stratified benchmark
(A), an adversarial benchmark (B), a three-method baseline comparison (C),
an identity×coverage threshold/sensitivity analysis (D), and the inputs
needed for real-data case studies (E). Experiment F (software engineering
evidence) is already complete and reported in `benchmarks/`, not here.

---

## 1. Data sources and access

| Resource | Endpoint | Version/date handling |
|---|---|---|
| UniProt | REST API, `reviewed:true` (Swiss-Prot only) | UniProt release recorded from the `X-UniProt-Release` response header at fetch time |
| RCSB PDB | sequence search API + polymer entity metadata | no fixed version; RCSB is continuously updated, so results are timestamped with fetch date, not a release number |
| PDBe SIFTS | `https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{accession}` | used as an independently-curated cross-reference source for Experiment E, separate from RCSB's own polymer-entity metadata |
| AlphaFold DB | prediction API (Section on retrieval in the manuscript) | queried live; a model's presence/absence is recorded as-observed at fetch time |
| ESM Atlas | fold API | queried live; sequences over `esm_max_length` (400 aa) are recorded as `invalid_input`, not attempted |

All fetches go through `bioseqinsight`'s own cache
(`~/.bioseqinsight/cache`), so a re-run within the cache TTL (7 days by
default) reproduces the same underlying API responses even if the live
databases have since changed. **Run `bioseqinsight cache clear` before the
main run if you want to guarantee fresh data**, and keep the default cache
otherwise if you want exact reproducibility of a specific run within the
TTL window.

## 2. Experiment A: stratified benchmark, inclusion criteria

**Target: 1,000 proteins** (500 is an acceptable fallback if time-constrained
— set `--target 500` and the same stratification proportions apply).

**Universe:** reviewed (Swiss-Prot) UniProt entries only. Unreviewed
(TrEMBL) entries are excluded because their annotation quality varies too
much for a ground-truth-adjacent benchmark.

**Length strata** (of the canonical UniProt sequence), each targeted at 25%
of the total:
- short: 50–150 aa
- medium: 150–300 aa
- long: 300–600 aa
- very long / complex: 600+ aa, OR annotated as having ≥2 domains (InterPro),
  OR annotated as having an isoform (UniProt `Alternative products` section)

**Within each stratum, deliberate sub-category sampling** (not simple random
sampling within the stratum), targeting roughly equal counts across the
sub-categories that are populated for that stratum:

1. **Exact/full-length** — has at least one PDB entry whose SEQRES covers
   ≥90% of the canonical sequence, OR has an AlphaFold DB model (used as the
   "exact" case for structures without full experimental coverage).
2. **Fragment** — has at least one PDB entry whose best chain covers <70% of
   the canonical sequence (construct/domain-only crystal structures).
3. **Missing-residue / partial** — has a PDB entry with unresolved residues
   recorded in the entry's own metadata (`pdbx_unobs_or_zero_occ_residues`
   or equivalent), i.e. the *right* protein with an incomplete model.
   **Currently unpopulated (see Amendment 1, below) — PDBe's
   `best_structures` endpoint does not carry per-residue gap information,
   and no verified alternative source has been wired in yet.**
4. **Isoform-bearing** — UniProt entry has ≥2 annotated isoforms.
5. **Paralogue family member** — UniProt entry belongs to a family with
   ≥2 annotated paralogues sharing >40% sequence identity (used later to
   construct Experiment B's within-family adversarial pairs; sampling here
   just ensures such families are represented).
6. **No structure anywhere** — no PDB entry, no AlphaFold model. Verified,
   not assumed: checked against all three resources at fetch time.

A protein may satisfy more than one category; category membership is
recorded per protein rather than forcing single-category assignment, and
sampling targets are treated as soft quotas, filled by random draw within
each (stratum × category) cell from the UniProt universe, seeded (Section 6).

**Exclusions:** membrane proteins with no structure at all in any resource
are capped at 10% of the "no structure" category (they are expected to
dominate it otherwise, which would bias that category away from soluble
proteins that genuinely lack coverage for other reasons). Sequences
containing non-standard amino acids beyond the 20 standard ones plus X are
excluded (parsing edge case, not a structure-retrieval question).

## 3. Experiment B: adversarial benchmark, construction

**Target: 200–300 cases.** Built from Experiment A's resolved candidates
plus deliberately constructed pairs, NOT sampled independently — every
adversarial case must have a documented, checkable reason for its expected
label.

| Class | Construction | Expected M-level (not enforced, evaluated) |
|---|---|---|
| B1 wrong paralogue | Query = protein A; structure forced = a structure whose owning accession is paralogue B (>40%, <95% identity to A) | M1 or M2, never M4 |
| B2 fragment | Query = full-length sequence; structure = a real PDB entry covering <70% of it | M2 or M3, never M4 (coverage-gated) |
| B3 wrong isoform | Query = isoform 1; structure = a structure whose SEQRES matches isoform 2 exactly | M1/M2 depending on isoform sequence divergence |
| B4 unrelated (negative control) | Query = random protein; structure = a randomly chosen, unrelated PDB entry | M0/M1 |
| B5 high-identity-wrong-entity | Query = protein A; structure = protein B with ≥95% identity to A but a different UniProt accession and different annotated function | Tests whether identity alone would wrongly justify M4; expected M3 at most (accession must also match for M4) |

Each case is a `(query_accession, query_sequence, forced_structure_source,
forced_identifier, expected_relationship)` tuple, stored in
`validation/data/adversarial_cases.csv`, generated once and frozen
(Section 6) before Experiment B is run.

## 4. Experiment C: baseline comparison, method definitions

Three methods, run over the **same** Experiment A candidate set, per
structural resource (RCSB and AlphaFold DB; ESM Atlas is excluded from this
comparison because it has no "native identifier retrieval" baseline — it
always predicts from sequence):

- **Method 1 — identifier-only:** for RCSB, take the first hit from a
  sequence search with no post-hoc filtering beyond the search's own
  ranking; for AlphaFold DB, this is definitionally equivalent to Method 2
  (accession lookup has no ranking step), and is reported as such rather
  than duplicated.
- **Method 2 — resource-native retrieval:** RCSB's own best-scoring hit
  after RCSB's internal identity/E-value filtering (i.e. what a user gets by
  trusting the search API's own default relevance ranking, one level more
  conservative than Method 1's "just take the top result"); AlphaFold DB by
  direct accession lookup.
- **Method 3 — `bioseqinsight` validated retrieval:** the full pipeline
  (Section 3.4–3.5 of the manuscript), M0–M4 assigned, only M4 counted as
  "accepted as correct."

**Ground truth for scoring Methods 1 and 2** (Section 5) is the same
adjudication used for Experiment A/E: independent SIFTS cross-reference
agreement, not `bioseqinsight`'s own alignment — see Section 5.

**Important asymmetry to report explicitly, not hide:** AlphaFold DB
retrieval is accession-keyed for all three methods, so Methods 1–3 cannot
differ in *what* is retrieved for AlphaFold DB, only in whether the M0–M4
layer additionally flags a low-coverage or low-identity result (which can
still happen — an AlphaFold model can be accession-correct but for a
different isoform than intended, for instance). Report RCSB and AlphaFold DB
comparisons **separately**; do not pool them into one false-association rate,
which would understate RCSB's error rate by diluting it with AlphaFold's
structurally-guaranteed accession matches.

## 5. Ground truth and avoiding circularity

Per the plan's Section 15, ground truth for "is this the correct structure"
is established independently of `bioseqinsight`'s own alignment engine,
using two signals that are both external to it:

1. **PDBe SIFTS residue-level UniProt mapping** — a curated, independently
   maintained cross-reference from PDB entry to UniProt accession and
   residue range, computed by EBI, not by this software.
2. **RCSB's own polymer-entity reference-sequence identifiers** — curated by
   RCSB at deposition, independent of any alignment `bioseqinsight` performs.

A case is labelled **ground-truth positive** (the structure IS the queried
protein) only if SIFTS and RCSB's own metadata agree that the entry's chain
maps to the queried accession. A case where these two independent sources
disagree with each other is labelled **ambiguous** and excluded from the
primary accuracy/precision/recall calculation, reported separately as its
own category (this is expected to be rare; if it exceeds 5% of cases,
investigate before proceeding, per the amendment rule above).

`bioseqinsight`'s own identity/coverage alignment is *not* part of ground
truth — it is the system under test. This is what Section 15 of the plan
calls avoiding circularity, and it is why SIFTS is fetched as a distinct
step (`validation/scripts/02_fetch_ground_truth.py`) before
`bioseqinsight` is ever run against a candidate.

## 6. Reproducibility

- Random seed: **20261001**, used for every stratified/random draw in this
  protocol. Changing it requires a logged amendment.
- The frozen accession list, once Experiment A's sampling completes, is
  written to `validation/data/dataset_accessions.csv` and is not
  re-sampled; a protein that later fails to resolve is recorded as
  unresolved, not silently replaced by a re-draw (silent replacement would
  let easy substitutions creep in).
- Every script records: `bioseqinsight` version, UniProt release header,
  fetch timestamps, and the exact command line used, into
  `validation/results/run_metadata.json`.
- API failures: a resource that fails after the configured retry budget
  (see `bioseqinsight`'s own fault-tolerance layer) is recorded as
  `service_error` for that candidate and resource, not retried indefinitely
  and not silently dropped from the denominator.

## 7. Metrics (Section 17 of the plan)

For every method (Experiment C) and the full system (Experiments A, B, D):
accuracy, precision, recall/sensitivity, specificity, F1, false-association
rate (= 1 − precision, restricted to cases where a structure was returned at
all), false-negative rate, Matthews correlation coefficient, and a full
confusion matrix, each with a 95% Wilson confidence interval for
proportions and a cluster (by query accession) bootstrap interval for
derived rates. Implemented in `validation/lib/stats.py`.

## 8. Pilot before the main run

Run `--target 100` first (Section 22 of the plan). Inspect
`validation/results/pilot_report.md` for: API errors, unexpectedly empty
categories, pathological sequences, and timing (to estimate the full run's
wall-clock time before committing to it). Do not change the protocol above
this line based on pilot *results* — only based on genuine mechanical
problems (e.g. a category that is empty because of a bug, not because it is
genuinely rare).

---

## Amendments log

**Amendment 1** (post-pilot, before any main run). The pilot run
(100 candidates) returned zero candidates in `exact_full_length`,
`fragment`, and `partial_missing_residue` despite ~29 candidates having a
UniProt PDB cross-reference. Cause: `lib/api.py`'s SIFTS integration called
`/mappings/uniprot/{accession}`, which is documented to take a **PDB ID** and
return its UniProt mapping (the reverse lookup), not a UniProt accession
returning its PDB structures; and its parsing additionally assumed a nested
`{accession: {"PDB": {pdb_id: [...]}}}` shape the real endpoint never
returns. This is a genuine mechanical defect (unverified assumption about
an external API's schema), not a disappointing result, and is corrected
rather than worked around: the correct, UniProt-accession-keyed endpoint is
PDBe's `best_structures` API (`/mappings/best_structures/{accession}`),
confirmed against two independent published client implementations,
returning a flat list per accession with a `coverage` field computed by
PDBe directly. Category 2 (fragment) and Category 1 (exact/full-length)
logic now uses this endpoint and schema. Category 3
(missing-residue/partial) is left unpopulated: no field in
`best_structures`'s response carries per-residue gap information, and
rather than repeat the same mistake by guessing a second unverified schema,
that category is documented as a known gap pending a verified data source.
A second, related bug was found and fixed at the same time:
`lib/ground_truth.py`'s `_sifts_agrees` (used by Experiment C/E, not
Experiment A) called the same wrong endpoint for a related but distinct
purpose (confirming whether one specific retrieved structure maps to one
specific accession); it now uses the PDB-ID-keyed
`/mappings/uniprot/{pdb_id}` endpoint, which is the documented, correct
tool for that specific question. This does not change Section 2's sampling
*targets* or *method*, only the functions supplying ground truth for three
of Experiment A's six categories and for Experiment C/E's adjudication, and
is logged here per this document's own rule that genuine mechanical fixes
are amended and re-run, not silently patched.

**Amendment 2** (post-pilot, after Amendment 1). Amendment 1's re-run pilot
confirmed the SIFTS fix (`exact_full_length`=37, `fragment`=22 in a
100-candidate pilot, both previously 0), but also showed `paralogue_family`
at 0. Cause: unlike the other five categories, `paralogue_family` was never
actually implemented in `classify_candidate()` -- Experiment B's adversarial
paralogue-pair construction (`01_build_adversarial.py`) has its own,
independent Pfam-based grouping and was unaffected, but Experiment A's
*sampling* had no signal to deliberately include paralogue-family members as
Section 2 Category 5 describes. Fixed by fetching `xref_pfam` as part of the
UniProt search result already made for every candidate (no extra network
calls) and adding a whole-pool tagging pass after classification: any
candidate sharing a Pfam ID with at least one other candidate in the pool is
tagged `paralogue_family` before stratified sampling runs. Logged here, not
silently applied, per this document's own rule.

**Amendment 3** (post-pilot, after Amendment 2). The pilot's 100 selected
proteins were all *Homo sapiens*, with accessions concentrated in the
alphabetically-first prefixes (O1: 21, O4: 13, O0: 11, P0: 10 ...). Section 2
specifies random sampling within strata, but `fetch_stratum_pool()` took the
first N results of a UniProt search, whose default order is not random, and
the seeded draw sampled only from that pool, so the seed controlled the
wrong step. Fixed by streaming every accession in each length range
(`/uniprotkb/stream`, `format=tsv`, verified live), drawing a seeded random
sample from the full list (seed derived from the global seed and the stratum
name), and fetching full entries for only those accessions. If the stream
fails, the old paging is used as a fallback and
`dataset_metadata.json` records `sampling_notes` saying the sample is biased,
plus `sampling_diagnostics` (distinct organisms, top-organism fraction, top
accession prefixes) for every run so concentration is always visible.
Consequence for reporting: the benchmark is a stratified, category-enriched
random sample of reviewed UniProt entries across all organisms, not a sample
of human proteins and not a representative sample of the proteome.

**Amendment 4** (after the first full run). Four mechanical defects, none a
disappointing result. (a) `rcsb_sequence_search` crashed on RCSB's HTTP 204
empty-body "no hits" reply (`json.loads("")`), aborting Experiment C at
protein 46 and leaving 306 of 782 adversarial cases unresolved; an empty body
now means no hits. (b) B2 (fragment) cases had no forced-structure resolution
path and all 68 were silently skipped; they now resolve to the query
accession's own PDB entry. (c) Paralogue pairs were banded on raw 4-mer
Jaccard 0.40-0.98 (~87-99.5% identity) instead of the intended 40-95%
identity, yielding 7 pairs against a target of 200-300; they are now banded
on estimated identity (J = p^4 / (2 - p^4), inverted) with a minimum
shared-k-mer count. (d) A single failing candidate aborted Experiment C; it
now becomes an error row. Not changed: the B3 (wrong isoform) label. In the
first run 17 of 79 scored B3 cases reached M4; this is not reported as a
validator failure, because an isoform differing from the canonical sequence
by a small region can yield a structure that legitimately matches the
canonical query at >=99% identity with the canonical accession
cross-referenced. B3 cases now record canonical and isoform lengths so their
label validity can be assessed before any B3 rate is reported. The B3 label
is a construction assumption, not independent ground truth.

**Amendment 5** (after the first full run). The first full run produced
results that must not be cited, because of defects in the analysis, not in the
data. (a) *Experiment C was a tautology.* Methods 1 and 2 were marked
"accepted" exactly when their hit was ground-truth positive and were then
scored against that same ground truth, giving precision = recall = 1.0 and
meaningless McNemar p-values. A naive method accepts whatever the resource
returns; correctness is now judged separately by independent ground truth,
the unit is the query protein, and proteins with a true structure (UniProt
lists a PDB entry) are reported apart from those without (where accepting
anything is a false association by definition). The validation layer is also
applied to the SAME naive hit, isolating its contribution. (b) *Experiment D
was circular.* Its "true positives" were bioseqinsight's own M4 decisions;
positives and negatives now come from independent ground truth and constructed
adversarial negatives, with an ablation showing what identity + coverage alone
would accept versus what the accession check removes. (c) *Ground truth checked
entity 1 only*; it now checks every polymer entity of the entry. (d)
*Experiment E asked for accessions,* which AlphaFold DB answers by
construction; it now adjudicates the top RCSB sequence-search hit. (e) *Finding
about the aligner:* unrelated pairs (B4) have median identity 28.6% at 84%
coverage, so the M2 identity threshold (30%) lies at the noise floor of this
gap-permissive alignment and 40 of 392 unrelated pairs (10.2%) reached M2. M2
therefore cannot be described as evidence of homology; the earlier suggestion
that those were chance homologues is withdrawn. The Methods statement that
percent identity is insensitive to the scoring scheme is withdrawn. M3/M4
(identity >= 95/99%) lie far above this floor.

**Amendment 6** (after inspecting the corrected run). Analysis-only. (a) The
shipped thresholds are M4 identity >= 99 / coverage >= 95, M3 95 / 80, M2
40 / 40; earlier documentation and the Experiment D script stated M3 95 / 90
and M2 30 / 50, which were wrong. (b) Experiment C now reports two operating
points, strict (M4) and permissive (M3 or better), and places every false
association in a tier by the hit's own identity/coverage (same sequence,
near-identical, different protein, unmeasured), because "false" means "not
this accession" and includes the same sequence filed under another accession.
(c) The cost of the strict gate is reported: why correct naive hits were not
accepted. (d) The unrelated-pair noise floor is reported conditioned on
coverage, because identity over the aligned region is inflated for short
subject chains. Corrected-run findings recorded here: naive top-hit retrieval
gave an accession-level false association for 730 of 1,000 queries (37 of
213 hits for proteins with a structure; 693 of 696 for proteins without);
validation-gated acceptance gave 0; the gate rejected 91 of 176 correct hits;
identity x coverage discriminated worse than identity alone because correct
hits include fragments; and all 7 v1.0-discordant queries now have a top RCSB
hit whose entry contains the queried accession.

**Amendment 7** (software, additive). Constructed fragments (B2) and 35 of 179
correct naive hits were confirmed, highly identical structures of the queried
protein that nonetheless fell to M1 or M2 because they cover only part of it,
and M1's message told the user the structure "does not correspond to the
requested protein". ValidationReport now carries `partial_structure` (accession
confirmed, identity >= 95, coverage < 80) and describe() reports it accurately.
The mapping level and the M4 gate are unchanged, so no reported
false-association or sensitivity figure changes. The module docstring no longer
asserts the v1.0 preproinsulin claim as fact. Known, unfixed limitation: the
RCSB provider downloads PDB-format files only, which do not exist for large or
new entries (23 of 179 correct naive hits and 123 of 730 naive false
associations could not be validated, all HTTP 404); validating from the entity
sequence when no PDB file exists is a bounded future fix.

**Amendment 8** (reproducibility and documentation). The adversarial set is now
frozen once used: `01_build_adversarial.py` refuses to overwrite an existing
file without `--force`, sorts the paralogue pairs before the seeded shuffle
(the earlier order depended on per-process string hashing, so the "seeded" draw
selected different pairs on different runs, and the set had been silently
rebuilt on each `run_all.sh adversarial`), and records a SHA-256. The file used
for the reported results was produced by the earlier code and is frozen as is.
README and docs now state the shipped thresholds (M3 95/80, M2 40/40). Every
number quoted in the manuscript's validation study is computed by
`07_manuscript_numbers.py` from the results files.

**Amendment 9** (analysis only). The first B3 validity rule (isoform more than 5%
shorter than the canonical sequence) tested the wrong quantity. Inspection of the
frozen results showed that all 7 such cases that reached M4 returned a chain of the
canonical length, i.e. the retrieval returned the full-length canonical structure,
so the "wrong isoform" label failed and M4 was correct; B3 is therefore not reported
as a failure rate and the rule is withdrawn as a validity criterion. The one B1 pair
that reached M4 has a chain of the same length as the query, 100% identical, with the
accession confirmed, so the forced structure is a structure of the query itself (a
construction error). The B2 fragments with a confirmed accession that are not flagged
partial all have measured identity below 95%. `07_manuscript_numbers.py` now computes
these diagnostics.
