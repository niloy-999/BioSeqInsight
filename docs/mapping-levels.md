# Sequence-to-structure mapping levels (M0–M4)

This is the reference for the classification that sits at the centre of
BioSeqInsight 2.0. If you read one document about this software, read this
one.

---

## The problem

You have a protein. You want its structure. You run a sequence search against
the PDB and something comes back. **Is it your protein?**

The question sounds trivial and is not. A sequence search returns things that
are *similar*, and similar covers:

* the same protein from a different species;
* a paralogue that shares 90 % of its sequence and none of its function;
* a fragment: one domain of a five-domain protein;
* a complex in which your protein is a minor component;
* a protein that binds or degrades yours, sharing local similarity in the
  interface;
* the same protein, correctly — which is what you wanted.

Software that reports "structure retrieved" for all six is not answering your
question. Version 1.0 of this software returned a hit for every query without being able
to say which case applied: for 7 of its 25 benchmark queries the accession
associated with the hit differed from the query. See
[`SOFTWARE_AUDIT.md`](SOFTWARE_AUDIT.md).

---

## The classification

Every structure BioSeqInsight retrieves is assigned one of five levels.

| Level | Name | Criteria | Interpretation |
|---|---|---|---|
| **M4** | Exact | Accession cross-referenced by the entry **and** identity ≥ 99 % **and** coverage ≥ 95 % | This is the structure of your protein. Cite it as such. |
| **M3** | High confidence | Identity ≥ 95 % and coverage ≥ 80 %, accession not confirmed | The same or a near-identical sequence, possibly filed under another accession. State the identity and coverage when you report it. |
| **M2** | Related record | Identity ≥ 40 % and coverage ≥ 40 %, below the M3 thresholds | A homologue, orthologue or fragment. Not a confirmed match for your protein (a confirmed but partial structure of it is flagged `partial_structure`, see below). Chance alone satisfies these thresholds for about 1 in 10 unrelated pairs, so do not read M2 as evidence of homology. |
| **M1** | Unverified | Either no comparison was possible, or a comparison was made and the structure is clearly not the query | No structural claim is supported. Read the message: the two cases are different. |
| **M0** | No result | No resource returned a structural record | — |

Thresholds are configurable (`exact_identity_threshold` and friends) and the
values in force are recorded with every result, so a reader can always see
which numbers produced a level.

---

## Why an accession match alone is not enough

An entry can cross-reference your accession and still not be a usable
structure of your protein: a 30-residue peptide from a 700-residue protein
legitimately carries the parent accession. Coverage catches this. That is why
M4 requires both the accession *and* the alignment, and why an accession
match with 45 % coverage lands in M2.

## Why a perfect sequence match alone is not enough

Two different proteins can share a near-identical region. More commonly, the
deposited construct is a chimera, a fusion, or an engineered variant whose
chain sequence matches your query over the aligned region while the entry
describes something else. Requiring the accession as well is what separates
"the sequences agree" from "the databases agree this is the same entity".

This is why a perfect 100 % identity hit with no accession confirmation is
M3, not M4. In practice M3 is usually correct. It is reported separately
because "usually correct" and "confirmed" are different statements, and the
difference matters when the result goes into a paper.

---

## How the measurements are made

**Alignment.** Semi-global (glocal) alignment: the query is aligned end to
end, the subject may be a sub-sequence without penalty. This is the right
model for comparing a full-length UniProt sequence with a deposited chain,
because a construct that omits a signal peptide and a disordered tail is
still an excellent structure of the folded part.

A strict global alignment would score that construct poorly and push a good
structure down to M2. A local alignment would score any shared domain highly
and push a bad one up to M4. Semi-global with a separate coverage term is the
combination that behaves correctly in both directions.

* **Identity** = matched positions / aligned positions × 100
* **Coverage** = aligned query positions / query length × 100

**Chain selection.** A coordinate file may hold several chains. The chain
that aligns best with the query is used; its identifier is reported so the
choice is visible.

**Accession cross-check.** Accessions are collected from the entry's
`DBREF` records, from `REMARK` cross-references, and from provider metadata
(for RCSB, the polymer entity's reference sequence identifiers). The query
accession is compared after normalisation, which handles FASTA header forms
like `sp|P24941|CDK2_HUMAN` and AlphaFold model names like `AF-P69905-F1`.

---

## Worked examples

**M4 — ubiquitin.** Query `P0CG48`. RCSB returns 1UBQ, whose polymer entity
cross-references the accession. Chain A aligns at 100 % identity over 100 %
coverage.

```
mapping        M4 Exact accession and sequence match
identity       100.0%   coverage 100.0%   chain A
```

**M2 — a bacterial homologue for a eukaryotic query.** The fold is shared,
the sequence is 40 % identical. Reported as a related record with a warning.
A reader who sees M2 and uses the structure anyway has made an informed
choice; one who was shown "structure retrieved" has not.

**M1 — a measured mismatch (format illustration).** When a comparison is made
and the best chain falls below the related thresholds, the report has this
shape. **The numbers are placeholders showing the format, not output from a
real query.**

```
mapping        M1 Unverified
WARNING: the retrieved structure does not correspond to the requested protein
(identity <X>%, coverage <Y>%), which is below the threshold for even a
related record. It is associated with <accession>.
```

Whether the historical preproinsulin query (`P01308`, associated with
`P14735` in the v1.0 benchmark) lands here, or in M4 because the returned
entry is an insulin-containing complex, is an empirical question answered by
`validation/scripts/06_case_studies.py`, not asserted here.

**M1 — a bare PDB identifier.** Querying `1UBQ` directly gives no query
sequence to compare against, so nothing can be verified. This is the "no
comparison possible" branch of M1, and the message says so. It is not a
failure; it is an honest absence of evidence.

**M0 — no coverage.** Bacteriophage T4 endolysin (`P00720`) had no AlphaFold
DB model at the time of the v1.0 benchmark. M0 with a note distinguishing a
coverage gap from a service failure — a distinction v1.0 did not draw.

---

## Reading a level in your own work

* **Writing a paper?** Only M4 supports "the structure of X (PDB 1ABC)".
  For M3, state the identity and coverage. For M2, describe it as a
  homologue and say which one.
* **Teaching?** M2 results are pedagogically useful precisely because they
  are wrong in an interesting way. Showing a class why a 40 %-identity hit is
  not their protein teaches more than a clean M4.
* **Screening many proteins?** The batch exporter includes the mapping level
  as a column. Filter on M4 before drawing any structural conclusion.

---

## Changing the thresholds

```bash
export BIOSEQINSIGHT_EXACT_IDENTITY_THRESHOLD=98
export BIOSEQINSIGHT_EXACT_COVERAGE_THRESHOLD=90
```

or in `~/.bioseqinsight/settings.json`. Settings are validated on load: the
related thresholds must not exceed the high-confidence ones, which must not
exceed the exact ones. The values used are written into every project and
every JSON export, so a lowered threshold cannot quietly inflate a result
table.

The defaults are deliberately strict. If a threshold has to be relaxed to
make a structure qualify as M4, that is information about the structure.


## The partial-structure qualifier

A deposited structure is often only a domain construct. If the structure's own
cross-references confirm the queried accession and the aligned region is at
least 95 % identical, but the structure covers less than 80 % of your sequence,
the result is flagged `partial_structure: true` and the message says so:

> Partial structure of the requested protein: the structure cross-references
> *accession* and matches the query at *x* % identity over the region it
> contains, but it covers only *y* % of your sequence. It is a real structure
> of this protein, not a full-length one.

The flag is a qualifier, not a level. The mapping level (M1 or M2) is unchanged
and the flag is never set for M3 or M4, so M4 still means a full-length exact
match.
