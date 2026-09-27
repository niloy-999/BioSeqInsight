# Tutorial 4 — Retrieving a structure and checking it is yours

**20 minutes. Network access required.**

This is the tutorial that covers what BioSeqInsight is for.

## 1. A clean case

```bash
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

Three things happened that you did not have to do:

1. The canonical sequence was fetched from UniProt, so there was something to
   compare against.
2. The deposited chain was aligned against it — 100 % identity, 100 %
   coverage.
3. The entry's UniProt cross-reference was checked against the accession you
   asked about.

Only when all three agree do you get **M4**. That is the one level that
supports writing "the structure of ubiquitin (PDB 1UBQ)" in a paper.

## 2. A case that is not clean

```bash
bioseqinsight structure P01308 --all
```

`P01308` is human preproinsulin. Watch what the sequence search returns, and
what BioSeqInsight says about it:

```
mapping        M1 Retrieved but unverified
WARNING: the retrieved structure does not correspond to the requested
protein (identity 8.2%, coverage 11.4%), which is below the threshold for
even a related record. It is associated with P14735.
```

`P14735` is insulin-*degrading* enzyme — the enzyme that destroys the
hormone. It is a legitimate sequence-search neighbour, because the enzyme's
structures contain bound insulin peptide. It is not insulin.

**Version 1.0 of this software reported this as a successful retrieval.** It
was one of 7 such cases in a 25-protein benchmark that claimed 25/25 success.
See [`../SOFTWARE_AUDIT.md`](../SOFTWARE_AUDIT.md).

## 3. The levels

| Level | Meaning | Safe to say "this is my protein"? |
|---|---|---|
| M4 | Accession confirmed **and** identity ≥ 99 % at ≥ 95 % coverage | **Yes** |
| M3 | Identity ≥ 95 %, coverage ≥ 90 %, accession unconfirmed | Probably — state the numbers |
| M2 | Homologue, orthologue or fragment | **No** |
| M1 | Not verified, or verified as different | **No** |
| M0 | Nothing returned | — |

Full detail: [`../mapping-levels.md`](../mapping-levels.md).

## 4. Watch the fallback work

```bash
bioseqinsight structure P0CG48 --all
```

`--all` queries every resource instead of stopping at the first exact match.
You will see RCSB, AlphaFold DB and ESM Atlas each reported separately, with
attempt counts:

```
  attempts       2 (#1 503, #2 200 ok)
  recovery       succeeded after an initial failure
```

When a service returns a transient error, the request is retried with
exponential backoff. When it stays down, the manager moves to the next
resource. When a resource simply has no model for your protein, that is
reported as a coverage gap — not as a failure.

## 5. Experimental versus predicted

RCSB gives experimental structures; the B-factor column holds
crystallographic temperature factors. AlphaFold DB and ESM Atlas give
predictions; the same column holds pLDDT confidence scores.

They occupy the same numeric range and mean completely different things.
BioSeqInsight labels them by provenance:

```
  mean pLDDT     92.41          # a predicted model
  mean B-factor  18.73          # an experimental structure
```

Version 1.0 guessed from the numeric range, which cannot distinguish a
well-ordered crystal structure from a confident prediction.

**Interpreting pLDDT:** above 90 very high, 70–90 confident, 50–70 low,
below 50 usually disordered rather than wrong. A low-pLDDT region is often a
correct statement that the region has no fixed structure.

## 6. Look at it

GUI → Structure tab → **View in 3Dmol.js** opens the structure in your
browser, with the mapping level shown alongside so you cannot forget whether
you are looking at the right molecule.

**Open external viewer** hands off to PyMOL or ChimeraX if installed.

## 7. Direct PDB access

```bash
bioseqinsight structure 1AQ1
```

```
mapping        M1 Retrieved but unverified
Structure retrieved; biological identity could not be verified.
```

You asked for an entry by name, so there was no query sequence to compare
against. That is not a failure — it is an honest absence of evidence, and it
is a different M1 message from the preproinsulin case. Read the message, not
just the level.

## 8. Cite the right people

BioSeqInsight retrieves structures; it does not predict them. If you use an
AlphaFold model, cite Jumper et al. (2021) and Varadi et al. (2022). For
ESMFold, cite Lin et al. (2023). For experimental structures, cite the entry
and Berman et al. (2000). The README has the full list.

## Next

[Tutorial 5 — Batch analysis and reproducible projects](05-batch-and-projects.md)
