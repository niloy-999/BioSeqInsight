# Example data

Small inputs for trying the software and for the tutorials in `docs/`.

| File | Use |
|---|---|
| `dna/demo_gene.fasta` | One clean ORF; the quickest way to see sequence analysis work |
| `dna/ambiguous.fasta` | Contains N codes at real positions; shows that ambiguity does not shift the reading frame |
| `dna/primer_panel.fasta` | Five short sequences spanning GC extremes, for melting-temperature comparison |
| `proteins/ubiquitin.fasta` | Resolves to an exact (M4) structural match |
| `proteins/benchmark_subset.fasta` | Five proteins covering an exact match, two near-identical paralogues, a low-confidence barrel, and the preproinsulin mis-mapping case |
| `proteins/disordered.fasta` | Intrinsically disordered; expect low prediction confidence |
| `structures/1UBQ.pdb` | An experimental coordinate file, for the local-file and viewer paths |

Try them:

```bash
bioseqinsight dna --file examples/dna/demo_gene.fasta
bioseqinsight protein --file examples/proteins/ubiquitin.fasta --sketch
bioseqinsight batch examples/proteins/benchmark_subset.fasta --formats csv,html
bioseqinsight structure P0CG48          # needs network access
```

`benchmark_subset.fasta` is deliberately chosen to produce different mapping
levels: run it with `--structures` and the exported table will show an M4
alongside cases that are not exact matches. That contrast is the point.
