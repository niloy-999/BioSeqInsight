# Tutorial 1 — Your first sequence analysis

**10 minutes. No network needed.**

## What you will do

Take a short gene, find its GC content and melting temperature, and
understand what the software is and is not telling you.

## 1. Run it

```bash
bioseqinsight dna --file examples/dna/demo_gene.fasta
```

```
Sequence: demo_gene
Length: 150 nt
GC content: 48.67%
Composition: A=41, T=36, G=39, C=34, N=0
Tm: 78.42 C (nearest_neighbor)
ORFs (>= 30 aa): 1
  strand + frame 1 nt 0-123 (41 aa)
```

## 2. Read the GC line carefully

```
Composition: A=41, T=36, G=39, C=34, N=0
```

The `N=0` is not decoration. Ambiguity codes are counted separately and
excluded from **both** the numerator and the denominator of the GC
calculation. A sequence that is half `N` would otherwise report an
artificially low GC content, because the unknown bases would be counted as
non-GC.

Try it:

```bash
bioseqinsight dna --file examples/dna/ambiguous.fasta
```

The `N=3` is reported, a warning appears, and the GC percentage describes
only the bases that are actually known.

## 3. Read the Tm line more carefully

```
Tm: 78.42 C (nearest_neighbor)
```

The method is stated because it matters. BioSeqInsight picks:

* **Wallace rule** for oligonucleotides up to about 20 nt — the familiar
  2(A+T) + 4(G+C);
* **SantaLucia (1998) nearest-neighbour** thermodynamics with salt
  correction for longer duplexes.

Compare a short primer with a long one:

```bash
bioseqinsight dna GCGCGGCCGCGATCGCGGCC
bioseqinsight dna ATTAATTAATATTATAAATT
```

The GC-rich primer melts far higher. That is real chemistry, not a quirk.

Now push past the valid range:

```bash
bioseqinsight dna --file examples/dna/demo_gene.fasta
```

A 150 nt sequence gets a nearest-neighbour value **and** a note that the
two-state model this assumes does not hold for long duplexes. Version 1.0 of
this software applied the Wallace rule to sequences of any length with no
comment at all, which produced large numbers with no physical meaning.

**The lesson:** a number without its method and conditions is not a result.
Everything BioSeqInsight reports carries both.

## 4. Get it as data

```bash
bioseqinsight dna --file examples/dna/demo_gene.fasta --json --out demo.json
```

```json
{
  "identifier": "demo_gene",
  "length": 150,
  "gc_percent": 48.67,
  "tm": {"tm_c": 78.42, "method": "nearest_neighbor", "na_mM": 50.0, "primer_nM": 250.0},
  "software_version": "2.0.0"
}
```

The conditions and the software version travel with the result. Six months
from now you will know exactly what produced that number.

## 5. Same thing in the GUI

```bash
bioseqinsight-gui
```

Sequence Analysis tab → open the FASTA file → **Full analysis**. Identical
values; the CLI and the GUI call the same functions.

## Next

[Tutorial 2 — ORFs, reading frames and ambiguity codes](02-orfs-and-translation.md)
