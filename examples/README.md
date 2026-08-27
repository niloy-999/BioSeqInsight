# Example sequences

## `short_dna.fasta`

Tiny teaching construct:

- Sequence: `CCCATGTCTTCTAAAGTTAAATAACCC`
- Frame +1 translation of the whole string starts at the leading `CCC` → `PMSSKVK` then `P`
- ORF search (ATG to stop) → `MSSKVK`
- Expected GC content: 37.04%
- Motif `ATG` 1-based start: 4

Use this on the **Sequence Analysis** tab.

## `ubiquitin.fasta`

Human ubiquitin monomer, 76 aa.

Approximate checks:

- Molecular weight ≈ 8565 Da
- Mean Kyte–Doolittle hydropathy ≈ −0.49

Use this for mass / hydropathy.

## `1UBQ.pdb`

Experimental ubiquitin from RCSB. Use **Load local PDB**. Crystal B-factors
are mobility values, not pLDDT.

New predictions and downloads are written to `structures/` as separate files.
