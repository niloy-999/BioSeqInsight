# Tutorial 6 — Using BioSeqInsight as a library

**20 minutes. Network optional.**

Every function returns data, not formatted text, so the package composes.

## 1. Sequence analysis

```python
from bioseqinsight import analyze_sequence

result = analyze_sequence(
    "ATGGCTAGCAAAGGTTTCCCGATTAAAGGCTTAACCGGATCAGCTTAA",
    identifier="demo",
)

print(result.gc_percent)          # 48.67
print(result.tm.tm_c, result.tm.method)
for orf in result.orfs:
    print(orf.strand, orf.frame, orf.aa_length, orf.protein)

import json
json.dumps(result.to_dict())      # everything is serialisable
```

## 2. Protein descriptors

```python
from bioseqinsight import ProteinAnalyzer

analyzer = ProteinAnalyzer(sequence, identifier="ubq")
result = analyzer.analyze(include_sketch=True)

print(result.molecular_weight, result.isoelectric_point, result.gravy)

# or individual functions
from bioseqinsight.core.protein import molecular_weight, hydropathy_profile
print(molecular_weight(sequence))
print(hydropathy_profile(sequence, window=9))
```

## 3. Structure retrieval with validation

```python
from bioseqinsight import StructureManager
from bioseqinsight.models.results import MappingLevel

outcome = StructureManager().retrieve("P24941")

if outcome.best_mapping is MappingLevel.M4_EXACT:
    best = outcome.best
    print(f"Confirmed: {best.identifier} at {best.download_path}")
else:
    print(f"Not confirmed: {outcome.best_mapping.label}")
    for warning in outcome.warnings:
        print(" ", warning)
```

**Always branch on the mapping level.** Treating any non-`None` result as
success is exactly the mistake this software exists to prevent.

## 4. Batch processing

```python
from bioseqinsight.workflows.batch import BatchRunner, BatchOptions
from bioseqinsight.io.export import write_csv

runner = BatchRunner()
records = runner.records_from_fasta_file("panel.fasta")
rows, summary = runner.run(
    records,
    BatchOptions(include_structures=True, workers=8),
    progress=lambda done, total, msg: print(f"{done}/{total}"),
)

print(summary.mapping_counts)
write_csv([r for r in rows if r.best_mapping == "M4"], "confirmed.csv")
```

## 5. Configure explicitly

```python
from bioseqinsight.config.settings import Settings
from bioseqinsight import StructureManager

settings = Settings.from_dict({
    "max_attempts": 5,
    "exact_identity_threshold": 98.0,
    "fallback_order": ("alphafold", "rcsb"),
    "download_dir": "./structures",
})
settings.validate()
manager = StructureManager(settings)
```

Settings are validated on load — inverted thresholds or an unknown provider
raise immediately rather than producing a subtly wrong result later.

## 6. Validate a structure you already have

```python
from pathlib import Path
from bioseqinsight.structures.validation import validate_structure
from bioseqinsight.config.settings import Settings

report = validate_structure(
    Path("1UBQ.pdb").read_text(),
    query_sequence=my_sequence,
    query_accession="P0CG48",
    settings=Settings(),
)
print(report.mapping_level.value, report.identity_percent, report.message)
```

Useful for auditing structures obtained some other way.

## 7. Test without a network

```python
from bioseqinsight.services.transport import FakeTransport, ScriptedResponse
from bioseqinsight.services.http_client import HttpClient, RetryPolicy
from bioseqinsight.structures.manager import StructureManager

transport = FakeTransport({
    "search.rcsb.org": [
        ScriptedResponse(503, "down"),          # first attempt fails
        ScriptedResponse(200, search_payload),  # retry succeeds
    ],
    "files.rcsb.org": [ScriptedResponse(200, coordinates)],
})
client = HttpClient(
    transport=transport,
    policy=RetryPolicy(max_attempts=3, backoff_base_s=0.0, jitter=False),
    sleeper=lambda s: None,
)
outcome = StructureManager(settings, client).retrieve("P0CG48")
assert outcome.best.recovered
```

This is how the whole test suite works, and how you can test code built on
BioSeqInsight without depending on a public service being up.

## 8. A complete script

```python
#!/usr/bin/env python3
"""Screen a FASTA file and report only confirmed structures."""
import sys
from bioseqinsight.workflows.batch import BatchRunner, BatchOptions
from bioseqinsight.io.export import write_csv

runner = BatchRunner()
rows, summary = runner.run(
    runner.records_from_fasta_file(sys.argv[1]),
    BatchOptions(include_structures=True, workers=6),
)

confirmed = [r for r in rows if r.best_mapping == "M4"]
print(f"{len(confirmed)}/{summary.total} confirmed at M4")
print("Mapping levels:", summary.mapping_counts)
write_csv(confirmed, "confirmed.csv")

for row in rows:
    if row.best_mapping in {"M1", "M2"}:
        print(f"  check {row.identifier}: {row.best_mapping}")
```

## API reference

Docstrings are the reference:

```python
help(analyze_sequence)
help(StructureManager.retrieve)
```

Every public function documents what it does **and** what it does not: which
formula, which valid range, and what the result must not be read as.
