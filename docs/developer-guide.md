# Developer guide

How to work on BioSeqInsight: the rules the codebase follows, how to extend
it, and what the test suite expects.

---

## Setting up

```bash
git clone https://github.com/niloy-999/BioSeqInsight.git
cd BioSeqInsight
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,validation,http]"
pre-commit install
python -m pytest -q
```

Installing `validation` and `http` locally is worth it: it means you run the
same tests the CI cross-validation and transport jobs run, rather than
discovering a difference after pushing.

---

## The rules

These are not style preferences. Each one exists because its absence caused a
specific problem in version 1.0.

**1. `core/` imports only the standard library.** No `requests`, no
Biopython, no numpy. A CI job installs the package bare and fails if any of
them can be imported. If you need a dependency for a calculation, the
calculation is in the wrong layer.

**2. Functions return data, not display strings.** Every public function
returns a value or a dataclass. Formatting belongs in `gui/formatting.py` or
`cli.py`. v1.0 returned formatted strings from every analysis function, which
made batch processing and export impossible without re-parsing text.

**3. Nothing in `core/` performs I/O.** No file reads, no network, no
`print`. `core` is pure computation, which is why it can be tested
exhaustively in milliseconds.

**4. No test touches the network.** Use `FakeTransport`. A test that depends
on RCSB being up is not a test, it is a monitor.

**5. Never substitute a different answer for a failed one.** If the requested
protein cannot be found, say so. Returning something similar and letting the
caller assume it is the same thing is the exact failure this rewrite exists
to prevent.

**6. Public functions have docstrings that say what the function does *and*
what it does not.** The `aliphatic_index` docstring says which formula;
`secondary_structure_sketch` says it is not a prediction. A CI job fails if a
module has no docstring.

---

## Repository layout

```
src/bioseqinsight/
├── core/         computation: alphabet, codon, translation, orf, motifs,
│                 thermodynamics, protein, sequence, alignment
├── models/       result dataclasses shared by every layer
├── services/     transport, http_client, cache, logging, errors
├── structures/   base, rcsb, alphafold, esmatlas, uniprot,
│                 pdbio, validation, manager, visualization
├── workflows/    batch, project
├── io/           export
├── gui/          formatting (Tkinter-free), widgets, views, main_window
└── cli.py
tests/            one module per source area, plus conftest fixtures
benchmarks/       data, scripts, protocols, results
docs/             this documentation
```

---

## Adding a structure provider

Providers are the most likely extension point. Each knows one resource and
nothing about the others.

**1. Subclass `StructureProvider`** in a new module under `structures/`:

```python
from .base import StructureProvider
from ..models.results import RetrievalStatus, StructureResult

class MyResourceProvider(StructureProvider):
    name = "myresource"

    def supports(self, query_type: str) -> bool:
        return query_type in {"uniprot", "sequence"}

    def retrieve(self, query, *, accession=None, sequence=None, refresh=False):
        if self.settings.offline:
            return self._skipped("offline mode")
        outcome = self.client.get(url, expect_text="ATOM")
        if not outcome.ok:
            return self._from_failed_request(outcome)
        path = self._save(outcome.text, identifier)
        return StructureResult(
            source=self.name,
            status=RetrievalStatus.OK,
            identifier=identifier,
            download_path=path,
            attempts=outcome.attempts,
            ...
        )
```

**2. Register it** in `manager.py`'s provider table and add its name to the
default `fallback_order` in `config/settings.py`.

**3. Do not validate inside the provider.** The manager calls
`validate_structure` on whatever you return. A provider that classified its
own results could exempt itself from the identity check, which defeats the
point.

**4. Write tests using `FakeTransport`.** Cover: success, 404, a transient
failure that recovers, a persistent failure, offline mode, and malformed
input. `tests/test_providers.py` has the pattern for each.

---

## Adding a calculation

Put it in the right `core/` module, return a value or a dataclass, and write
tests that cover:

* at least one value verifiable against an independent published source;
* the empty input;
* input containing ambiguity codes;
* the boundary where the method stops being valid, with a warning attached.

If Biopython implements the same thing, add a parity case to
`tests/test_biopython_parity.py`. If your implementation deliberately differs
from Biopython's, the test documents the difference and why — see
`TestDocumentedDifferences` for the Tm and ambiguity cases.

---

## Testing

```bash
python -m pytest -q                        # whole suite, offline, ~20 s
python -m pytest tests/test_validation.py -v
python -m pytest --cov=bioseqinsight --cov-report=term-missing
python -m pytest -k "mapping or identity"
```

**Fixtures** (in `conftest.py`): `tmp_settings` confines every path to a
temporary directory and disables backoff jitter; `working_transport` scripts
a healthy response from every bundled service; `ubiquitin_pdb` is a real
coordinate file; `make_client` wires an `HttpClient` to a fake transport with
instant backoff.

**Writing a failure-path test:**

```python
def test_recovery_after_transient_failure(self, tmp_settings):
    transport = FakeTransport({"esmatlas.com": [
        ScriptedResponse(504, "gateway timeout"),
        ScriptedResponse(200, coordinates),
    ]})
    result = provider(transport, tmp_settings).retrieve(SEQ, sequence=SEQ)
    assert result.status is RetrievalStatus.OK
    assert result.recovered is True
```

The last scripted response repeats once the list is exhausted, so a
single-element list models a persistently failing service.

**Coverage.** The four Tkinter view modules are omitted from coverage in
`pyproject.toml` because they cannot run on a headless runner. This is why
the formatting logic they display lives in `gui/formatting.py`, which is
covered at 95 % — the wording a user reads is tested even though the widgets
are not.

---

## Style and static analysis

```bash
ruff check src tests benchmarks
black src tests benchmarks
mypy src/bioseqinsight
```

Line length 100. `pre-commit install` runs all of this plus the test suite
before each commit.

`mypy` is advisory in CI while type coverage is completed. New code should be
annotated.

---

## Releasing

1. Update `CHANGELOG.md`.
2. Bump the version in `src/bioseqinsight/__init__.py`, `pyproject.toml` and
   `CITATION.cff`. A CI job fails if the three disagree.
3. `python -m pytest --cov=bioseqinsight --cov-fail-under=85`
4. Re-run the offline benchmarks and commit the results.
5. Re-run the live structure benchmark and commit its dated output.
6. Tag `vX.Y.Z` and push. The release workflow verifies on three platforms,
   builds, checks the wheel and drafts a GitHub release.
7. Publish the release; the Zenodo integration mints the DOI.

---

## Performance notes

Local analysis runs at roughly 850–900 sequences per second on a modern
x86-64 machine, and per-sequence cost is flat from 10 to 1 000 sequences.

The cost centre is `core/alignment.py`: Needleman–Wunsch is O(mn) in time and
memory, which is why there is a 4-million-cell cap that raises
`AlignmentTooLarge` rather than exhausting memory. Two 3 000-residue
sequences exceed it deliberately. If you need to align longer sequences, the
right answer is a different algorithm, not a bigger cap.

Batch structure retrieval is I/O-bound and uses a thread pool. Local-only
batches are CPU-bound and single-threaded, which is why the published
throughput figure is honest rather than flattering.
