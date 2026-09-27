# Tutorial 5 — Batch analysis and reproducible projects

**20 minutes. Network optional.**

## 1. Run a file

```bash
bioseqinsight batch examples/proteins/benchmark_subset.fasta --formats csv,html
```

Five proteins, one table. The summary reports how many succeeded, how long it
took and the throughput.

Open the HTML report in a browser. It is self-contained — no external CSS, no
external scripts — so it can be emailed, archived, or opened in five years
and still render.

## 2. Failures do not stop the run

Make a file with one broken record:

```bash
printf '>good\nATGAAACCCGGGTAA\n>broken\n@@@@@@@@\n>alsogood\nATGGCTTAA\n' > mixed.fasta
bioseqinsight batch mixed.fasta --formats csv
```

Three rows out. The broken one carries its error in the `error` column; the
other two are analysed normally. A thousand-sequence job completes and tells
you which three inputs were malformed, rather than dying on record 417.

## 3. Add structures

```bash
bioseqinsight batch examples/proteins/benchmark_subset.fasta \
    --structures --workers 4 --formats csv
```

Much slower, and it needs network access. The summary now includes structures
attempted and retrieved, cache hits, recovered requests, and the distribution
of mapping levels.

**Filter on the mapping level before drawing any conclusion:**

```bash
awk -F, 'NR==1 || $0 ~ /,M4,/' bioseqinsight_batch.csv > confirmed.csv
```

`benchmark_subset.fasta` is chosen so this matters: it contains an exact
match, two near-identical paralogues (HBA and HBB), a beta-barrel, and the
preproinsulin case. They will not all come back M4, and that contrast is the
point of the file.

## 4. Projects

```bash
bioseqinsight project create ./study --name "Kinase panel" \
    --fasta examples/proteins/benchmark_subset.fasta
bioseqinsight batch examples/proteins/benchmark_subset.fasta \
    --project ./study --structures
bioseqinsight project info ./study
```

```
study/
├── project.json      manifest, record count, notes
├── metadata.json     software version + complete settings snapshot
├── sequences/        the exact inputs analysed
├── structures/       retrieved coordinates
├── results/          exported tables
└── logs/             structured run logs
```

Look at `metadata.json`. It holds the software version and **every** setting
in force: retry budget, mapping thresholds, Tm conditions, genetic code,
endpoints. Six months from now, "what thresholds did I use?" has an answer.

## 5. Hand it over

```bash
bioseqinsight project export ./study --out study.zip
bioseqinsight project import study.zip ./restored
bioseqinsight project info ./restored
```

The archive re-imports into a working project. That is the difference between
sending a colleague a spreadsheet and sending them something they can check.

Extraction refuses entries whose paths would escape the destination
directory, so importing an archive from someone else cannot write outside the
folder you chose.

## 6. Notes

```python
from bioseqinsight.workflows.project import Project
project = Project.open("./study")
project.write_note("Re-ran after UniProt release 2026_03; P04637 changed length")
```

Timestamped, stored in the manifest, travels with the export.

## Next

[Tutorial 6 — Using BioSeqInsight as a library](06-scripting.md)
