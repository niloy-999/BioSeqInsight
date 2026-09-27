# User guide

Covers both interfaces. They sit on the same core API, so anything here
applies to whichever you use.

---

## Contents

1. [Graphical interface](#graphical-interface)
2. [Command-line interface](#command-line-interface)
3. [Batch analysis](#batch-analysis)
4. [Projects](#projects)
5. [Caching](#caching)
6. [Reading structure results](#reading-structure-results)
7. [Common tasks](#common-tasks)

---

## Graphical interface

```bash
bioseqinsight-gui
```

Four tabs, a status bar, and menus for projects, tools and help.

### Sequence Analysis

Paste a nucleotide sequence or open a FASTA file. FASTA headers, blank lines,
position numbers and whitespace are handled; you do not need to clean the
input first.

| Button | Result |
|---|---|
| **Full analysis** | Length, GC content, composition, melting temperature, ORFs and, if a motif is entered, its hits |
| **GC content** | GC percentage with the counts behind it; ambiguity codes are excluded from both terms |
| **Reverse complement** | The reverse complement, wrapped |
| **Transcribe to mRNA** | T → U |
| **Six-frame translation** | All six frames, `*` marking stops |
| **Find ORFs** | ORFs above the minimum length, both strands |
| **Melting temperature** | Tm with the method and conditions used |
| **Find motif** | IUPAC motif search on both strands |

**Minimum ORF length** defaults to 30 aa. Lower it for short constructs.
**Motif** accepts IUPAC degenerate codes: `GAATTC` is exact, `WGATAR` matches
`[AT]GATA[AG]`.

Export the current result as CSV or JSON with the buttons at the bottom.

### Protein Analysis

Molecular weight, GRAVY, isoelectric point, aromaticity, aliphatic index,
extinction coefficients (reduced and with cystines) and composition.

**Hydropathy profile** draws a Kyte–Doolittle profile as a text chart; the
window size is adjustable and must be odd.

**Propensity sketch** shows helix/sheet/coil propensities. Read the label:
this is an illustrative sketch, not a structure prediction. It is included
because it is pedagogically useful and labelled because it is easy to
misread.

**Send to Structure tab** copies the sequence across, which is the usual next
step.

### Structure Retrieval

Accepts three kinds of input, detected automatically:

* a protein sequence (searched by sequence, or folded by ESM Atlas);
* a UniProt accession such as `P24941` (the canonical sequence is fetched
  first, so the identity check has something to compare against);
* a PDB identifier such as `1AQ1` (downloaded directly).

| Button | Behaviour |
|---|---|
| **Retrieve and validate** | Try resources in order, stop at an exact (M4) match |
| **Query all resources** | Query every resource regardless, for comparison |
| **AlphaFold DB / RCSB / ESM Atlas only** | Restrict to one resource |
| **Load local PDB** | Inspect a file you already have |
| **View in 3Dmol.js** | Open the structure in your browser |
| **Open external viewer** | Hand off to PyMOL, ChimeraX or similar |

The coloured badge above the result is the mapping level. The output lists
every resource tried, how many attempts each took, whether it recovered from
a failure, whether it came from cache, and the identity and coverage behind
the verdict.

The window stays responsive during retrieval: network calls run on a
background thread.

### Batch Analysis

Choose a FASTA file, set the options, click **Run batch**. The log reports
progress and a summary; export to CSV, TSV, JSON or HTML.

Tick **Retrieve and validate structures** to include structure lookup. This
is much slower and needs network access; workers default to 4.

---

## Command-line interface

```bash
bioseqinsight --help
bioseqinsight <command> --help
```

Global options: `--offline`, `--no-cache`, `--refresh`, `--download-dir`,
`--config`, `--log-level`, `--quiet`.

Exit codes: `0` success, `2` user error, `3` runtime failure.

### dna

```bash
bioseqinsight dna ATGGCTAGC...
bioseqinsight dna --file gene.fasta --min-orf 50 --motif GAATTC
bioseqinsight dna --file gene.fasta --json --out result.json
cat gene.fasta | bioseqinsight dna -
```

`--table` selects the genetic code (1 standard, 2 vertebrate mitochondrial,
4 mould/protozoan mitochondrial, 11 bacterial).

### protein

```bash
bioseqinsight protein MQIFVKTLTGK... --sketch
bioseqinsight protein --file protein.fasta --json
```

### structure

```bash
bioseqinsight structure P24941
bioseqinsight structure 1AQ1
bioseqinsight structure MQIFVKTLTGK... --providers rcsb,alphafold
bioseqinsight structure P24941 --all --json --out cdk2.json
```

Exits `3` when nothing was retrieved, so it composes in scripts:

```bash
if bioseqinsight structure "$ACC" --json --out "$ACC.json"; then
    echo "found"
fi
```

### batch

```bash
bioseqinsight batch sequences.fasta
bioseqinsight batch sequences.fasta --structures --workers 8 \
    --formats csv,json,html --out-prefix results/run01
bioseqinsight batch sequences.fasta --project ./my_project
```

### project, cache, config, version

```bash
bioseqinsight project create ./study --name "Kinase panel" --fasta kinases.fasta
bioseqinsight project info ./study
bioseqinsight project export ./study --out study.zip
bioseqinsight project import study.zip ./restored

bioseqinsight cache stats
bioseqinsight cache list --namespace rcsb
bioseqinsight cache clear

bioseqinsight config
bioseqinsight config --save ~/.bioseqinsight/settings.json
bioseqinsight version
```

---

## Batch analysis

One record failing does not stop the run. Failures appear as rows with an
`error` column, so a thousand-sequence job completes and tells you which
three inputs were malformed.

The CSV columns are fixed and documented in `io/export.py` (`BATCH_COLUMNS`),
so downstream scripts can rely on them. The HTML report is self-contained —
no external CSS or scripts — which means it can be emailed or archived and
will still render years later.

---

## Projects

A project makes a run repeatable:

```
study/
├── project.json      manifest and notes
├── metadata.json     software version + complete settings snapshot
├── sequences/        the exact inputs
├── structures/       retrieved coordinates
├── results/          exported tables
└── logs/             structured logs
```

```bash
bioseqinsight project create ./study --fasta input.fasta
bioseqinsight batch input.fasta --project ./study --structures
bioseqinsight project export ./study --out study.zip
```

The exported zip re-imports into a working project, which is what makes
"here is everything needed to check my analysis" one action.

---

## Caching

Responses are cached on disk with a SHA-256 digest and a TTL (7 days by
default). A tampered or corrupt entry is discarded rather than served.

```bash
bioseqinsight cache stats           # hit rate, entry count, size
bioseqinsight --refresh structure P24941    # ignore the cache this once
bioseqinsight --no-cache structure P24941   # bypass entirely
```

Caching is why re-running a benchmark is fast and why it does not hammer
public services.

---

## Reading structure results

```
[RCSB] ok
  identifier     1UBQ
  mapping        M4 Exact accession and sequence match
  identity       100.0%   coverage 100.0%   chain A
  attempts       2 (#1 503, #2 200 ok)
  recovery       succeeded after an initial failure
```

* **mapping** is the verdict. See [`mapping-levels.md`](mapping-levels.md).
* **identity / coverage** are the measurements behind it.
* **attempts** shows the retry history; **recovery** appears when a request
  succeeded only after failing.
* **cache** appears when the payload came from disk.

Only **M4** supports stating that the structure is your protein.

---

## Common tasks

**Characterise a gene and get its product's structure**

```bash
bioseqinsight dna --file gene.fasta --json --out gene.json
python -c "import json;print(json.load(open('gene.json'))['orfs'][0]['protein'])" > peptide.txt
bioseqinsight structure "$(cat peptide.txt)"
```

**Screen a panel and keep only confirmed structures**

```bash
bioseqinsight batch panel.fasta --structures --formats csv
awk -F, 'NR==1 || $0 ~ /,M4,/' bioseqinsight_batch.csv > confirmed.csv
```

**Compare melting temperatures across primers**

```bash
bioseqinsight batch examples/dna/primer_panel.fasta --formats csv
```

**Work entirely offline**

```bash
bioseqinsight --offline batch sequences.fasta
```

Local analysis runs normally; structure retrieval is skipped and reported as
skipped rather than failing silently.
