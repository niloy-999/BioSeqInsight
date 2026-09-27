# Tutorial 3 — Protein properties and what they mean

**15 minutes. No network needed.**

## 1. Analyse ubiquitin

```bash
bioseqinsight protein --file examples/proteins/ubiquitin.fasta --sketch
```

```
Length: 76 aa
Molecular weight: 8564.74 Da (8.565 kDa)
GRAVY: -0.4895
Isoelectric point: 6.79
Aromaticity: 0.0395
Aliphatic index: 100.00
Extinction (280 nm): 1490 reduced / 1490 with cystines
```

The accepted average mass of ubiquitin is 8564.8 Da. Agreement to 0.06 Da is
a useful sanity check that the residue masses and water subtraction are
right — and it is asserted as a test
(`test_ubiquitin_matches_the_published_value`), so a regression would break
the build rather than reach you.

## 2. What each number means

**Molecular weight** — sum of average residue masses minus one water per
peptide bond. Use it to predict gel migration or check mass spectrometry.

**GRAVY** — mean Kyte–Doolittle hydropathy. Negative is hydrophilic, positive
hydrophobic. Ubiquitin's −0.49 fits a soluble cytosolic protein; a value
above about +0.4 suggests membrane association.

**Isoelectric point** — the pH where net charge is zero. Sequence-based pI
ignores local structural effects and buried residues, as every sequence-based
calculator does. It is a starting point for choosing a buffer, not a
measurement.

**Extinction coefficient** — Gill & von Hippel, from Trp, Tyr and cystine
content. Use it to convert A280 to concentration. The two values differ when
cysteines form disulfide bonds; ubiquitin has no cysteines, so they match.

**Aliphatic index** — relative volume in aliphatic side chains; correlates
with thermostability. Ubiquitin's 100 is high, consistent with a notably
stable protein.

## 3. Hydropathy profile

In the GUI: Protein tab → **Hydropathy profile**. A window of 9 is standard
for general hydrophobicity; 19 approximates a membrane-spanning helix.

Compare a soluble protein with a disordered one:

```bash
bioseqinsight protein --file examples/proteins/disordered.fasta
```

Alpha-synuclein has a strongly negative GRAVY and a biased composition — the
sequence-level signature of intrinsic disorder.

## 4. The propensity sketch, and what it is not

```bash
bioseqinsight protein --file examples/proteins/ubiquitin.fasta --sketch
```

```
Secondary-structure propensity sketch
  helix 34.2%   sheet 28.9%   coil 36.8%
  NOT a structure prediction; do not report these percentages as such.
```

This applies Chou–Fasman-style per-residue propensities. It is a 1970s
statistical method, useful for teaching and for a rough feel, and it is
**not** a structure prediction. Modern predictors use evolutionary profiles
and neural networks and are far more accurate.

The disclaimer is printed every time, in the GUI, the CLI and the exported
data, because version 1.0's wording invited exactly this misreading. If you
need real secondary structure: retrieve a structure (Tutorial 4) and run
DSSP, or use a dedicated predictor.

## 5. Cross-validation

If you have Biopython installed:

```bash
pip install biopython
python -m pytest tests/test_biopython_parity.py -v
```

Molecular weight, GRAVY, aromaticity, isoelectric point, translation, GC
content and reverse complement are all compared against Biopython. Where the
two deliberately differ — melting-temperature salt correction defaults,
ambiguity handling — the test documents why rather than asserting equality.

BioSeqInsight implements these itself so that it has no required
dependencies. The parity suite is what turns that from a claim into a
demonstration.

## Next

[Tutorial 4 — Retrieving a structure and checking it is yours](04-structure-retrieval.md)
