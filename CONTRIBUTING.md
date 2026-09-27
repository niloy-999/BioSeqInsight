# Contributing to BioSeqInsight

Contributions are welcome — bug reports, corrections to the science,
documentation, and code.

## Before you start

Read [`docs/architecture.md`](docs/architecture.md) and the "rules" section
of [`docs/developer-guide.md`](docs/developer-guide.md). They are short, and
each rule exists because its absence caused a specific problem in version
1.0.

## Setting up

```bash
git clone https://github.com/niloy-999/BioSeqInsight.git
cd BioSeqInsight
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,validation,http]"
pre-commit install
python -m pytest -q
```

## Reporting a bug

Include the output of `bioseqinsight version`, the exact command or code, the
input (a minimal FASTA record if possible), what you expected and what
happened.

**Scientific disagreements are bug reports and are especially welcome.** If a
calculation disagrees with a published value or with another tool, say which
value you expected and where it comes from. Several of the defects fixed in
2.0 were of exactly this kind.

## Requesting a feature

Say what you are trying to do, not just what to add. Requests that would make
the software imply more certainty than it has — for example, removing the
"not a prediction" label from the propensity sketch, or auto-accepting M2
matches — will be declined, with an explanation.

## Pull requests

1. Branch from `develop`.
2. Write the test first where you can. New behaviour needs a test; a bug fix
   needs a test that fails before the fix.
3. Keep the suite offline and deterministic. Use `FakeTransport`; never call
   a real service from a test.
4. Run `python -m pytest -q`, `ruff check src tests benchmarks`,
   `black src tests benchmarks`.
5. Update `CHANGELOG.md` under "Unreleased".
6. Update the docs if behaviour changed.

Small, focused pull requests get reviewed faster than large ones.

## What a good test looks like

Tests here assert on values and on the wording users see, not on
implementation details:

```python
def test_ubiquitin_matches_the_published_value(self):
    # The accepted average mass of ubiquitin is 8564.8 Da.
    assert molecular_weight(UBIQUITIN) == pytest.approx(8564.8, abs=0.5)
```

The comment cites where the expected value comes from. Do that.

For anything touching structure retrieval, cover the failure paths: a
transient error that recovers, a persistent one, a 404, and offline mode.

## Adding a calculation

Cover the normal case against an independently published value, the empty
input, ambiguity codes, and the boundary where the method stops being valid
(with a warning attached to the result). If Biopython implements the same
thing, add a case to `tests/test_biopython_parity.py`.

## Adding a structure provider

See the walkthrough in [`docs/developer-guide.md`](docs/developer-guide.md).
The one rule that is not negotiable: providers do not validate their own
results. The manager calls `validate_structure` on whatever a provider
returns, so no resource can exempt itself from the identity check.

## Documentation

Prose, not bullet fragments. Say what something does and what it does not.
If you document a number, say which script produced it.

## Code of conduct

Be civil and assume good faith. Disagree about the science as much as you
like; do not be unpleasant about it. Report problems through GitHub's private
reporting to the maintainers.

## Licence

Contributions are accepted under the MIT licence of the project.
