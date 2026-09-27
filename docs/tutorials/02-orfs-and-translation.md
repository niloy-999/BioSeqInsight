# Tutorial 2 — ORFs, reading frames and ambiguity codes

**15 minutes. No network needed.**

## What you will do

Find open reading frames, and see why one apparently harmless shortcut —
deleting ambiguity codes before translating — produces confidently wrong
answers.

## 1. Find the ORFs

```bash
bioseqinsight dna --file examples/dna/demo_gene.fasta --min-orf 30
```

Each ORF reports its strand, frame, nucleotide coordinates and length. The
minimum length filters noise: in random sequence, short ORFs appear by
chance, and 30 aa is a common threshold for a real coding region.

Lower it and watch spurious ORFs appear:

```bash
bioseqinsight dna --file examples/dna/demo_gene.fasta --min-orf 5
```

Coordinates are 0-based and half-open, and reverse-strand ORFs are reported
in forward-strand coordinates so they can be compared directly.

## 2. Six frames at once

```bash
bioseqinsight-gui       # Sequence tab -> Six-frame translation
```

Three forward frames, three reverse, with `*` marking stops. The frame with
long stretches between stops is usually the coding one.

## 3. The ambiguity trap

Real sequencing data contains `N`. What should translation do with it?

```bash
bioseqinsight dna --file examples/dna/ambiguous.fasta --min-orf 5 --json \
  | python -c "import json,sys; print(json.load(sys.stdin)['orfs'][0]['protein'])"
```

```
MKXGFPIKGLTGSAYDPKGLADLKAWIERLGKPVKAA
```

Note the `X` at position 3. The codon containing `N` is unknown, so it
translates to `X`, **and the codon stays in place**.

Here is why that matters. If ambiguity codes were deleted first:

```
ATGAAANNNGGTTTC...      12 bases consumed as 4 codons -> M K X G
ATGAAAGGTTTC...         after deleting NNN: M K G F ...
```

Deleting three bases happens to keep the frame here. Delete **one**:

```bash
python -c "
import sys; sys.path.insert(0, 'src')
from bioseqinsight.core.translation import translate
print(translate('ATGNAAAGGTTT'))
"
```

```
MXRF
```

Correct. Had the `N` been stripped, `ATGAAAGGTTT` would translate to `MKG` —
a plausible-looking peptide that is wrong from residue two onward, with
nothing to indicate a problem.

**This was a real defect in version 1.0**, found by writing the test that is
now `test_ambiguity_codes_translate_to_x`. It is the clearest example in this
codebase of a bug that produces believable output instead of an error, which
is the kind that survives.

## 4. Genetic code tables

```bash
bioseqinsight dna ATGAAACCCGGGTGATAA --table 1 --min-orf 2   # TGA = stop
bioseqinsight dna ATGAAACCCGGGTGATAA --table 2 --min-orf 2   # TGA = Trp
```

Table 1 is standard, 2 is vertebrate mitochondrial, 4 mould/protozoan
mitochondrial, 11 bacterial. Using the wrong table on mitochondrial or
bacterial sequence gives premature stops and truncated ORFs.

## 5. Motifs

```bash
bioseqinsight dna --file examples/dna/demo_gene.fasta --motif GAATTC
bioseqinsight dna --file examples/dna/demo_gene.fasta --motif WGATAR
```

`WGATAR` expands to `[AT]GATA[AG]` — full IUPAC degeneracy. Both strands are
searched and hits are reported in forward coordinates.

This is exact pattern matching, not statistical motif discovery. For
discovering unknown motifs, use MEME.

## Next

[Tutorial 3 — Protein properties and what they mean](03-protein-properties.md)
