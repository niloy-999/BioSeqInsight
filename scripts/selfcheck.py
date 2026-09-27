#!/usr/bin/env python3
"""Post-installation self-check.

Verifies that an installed copy of BioSeqInsight is sound: that the package
imports, that the calculations agree with independently published reference
values, that the retry machinery behaves, and that file output works. It
needs no network and no test framework, so it can be run on a machine where
the repository is not available.

    python scripts/selfcheck.py
    python -m bioseqinsight.selfcheck     # if installed from a checkout

Exit code 0 means every check passed.
"""

from __future__ import annotations

import json
import sys
import tempfile
import traceback
from pathlib import Path

try:  # running from a source checkout
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
except Exception:  # pragma: no cover
    pass

PASSED: list[str] = []
FAILED: list[tuple[str, str]] = []

UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)


def check(name: str):
    """Decorator registering a check and recording its outcome."""

    def wrap(function):
        try:
            function()
        except Exception as exc:
            FAILED.append((name, f"{type(exc).__name__}: {exc}"))
            print(f"  FAIL  {name}")
            if "-v" in sys.argv:
                traceback.print_exc()
        else:
            PASSED.append(name)
            print(f"  ok    {name}")
        return function

    return wrap


print("BioSeqInsight self-check\n" + "=" * 60)

# -- 1. the package imports -------------------------------------------------
print("\nInstallation")


@check("package imports")
def _import():
    import bioseqinsight

    assert bioseqinsight.__version__


@check("core has no third-party dependencies")
def _no_deps():
    import subprocess

    # Block the optional third-party modules, then import the core. If any
    # core module has grown an accidental dependency, this fails.
    paths = json.dumps([p for p in sys.path if p])
    code = (
        f"import sys; sys.path[:0] = {paths};"
        "sys.modules['requests'] = None; sys.modules['Bio'] = None;"
        "import bioseqinsight.core.sequence, bioseqinsight.core.protein,"
        " bioseqinsight.core.alignment, bioseqinsight.core.thermodynamics;"
        "print('ok')"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr.strip()


# -- 2. the science ---------------------------------------------------------
print("\nCalculations (against independently published values)")


@check("ubiquitin molecular weight is 8564.8 Da")
def _mw():
    from bioseqinsight.core.protein import molecular_weight

    value = molecular_weight(UBIQUITIN)
    assert abs(value - 8564.8) < 0.5, value


@check("GC content of GGCC is 100%")
def _gc():
    from bioseqinsight.core.sequence import gc_content

    assert gc_content("GGCC") == 100.0
    assert gc_content("AATT") == 0.0


@check("reverse complement round-trips")
def _rc():
    from bioseqinsight.core.codon import reverse_complement

    sequence = "ATGCGATCGTAGCTAGCTA"
    assert reverse_complement(reverse_complement(sequence)) == sequence
    assert reverse_complement("GAATTC") == "GAATTC"


@check("standard genetic code has 64 codons and 3 stops")
def _code():
    from bioseqinsight.core.codon import get_genetic_code

    code = get_genetic_code(1)
    assert len(code.forward) == 64
    assert code.stops == frozenset({"TAA", "TAG", "TGA"})


@check("ambiguity codes do not shift the reading frame")
def _frame():
    from bioseqinsight.core.translation import translate

    # The v1.0 defect: stripping N shifted every downstream codon.
    assert translate("ATGNAAAGGTTT") == "MXRF"


@check("mitochondrial code differs from standard where it should")
def _mito():
    from bioseqinsight.core.translation import translate

    assert translate("TGA", table=1) == "*"
    assert translate("TGA", table=2) == "W"


@check("GC-rich duplex melts higher than AT-rich")
def _tm():
    from bioseqinsight.core.thermodynamics import nearest_neighbor_tm

    gc = nearest_neighbor_tm("GCGCGCGCGCGCGCGCGCGC")
    at = nearest_neighbor_tm("ATATATATATATATATATAT")
    assert gc > at + 20, (gc, at)


@check("ORF detection finds a known reading frame")
def _orf():
    from bioseqinsight.core.orf import find_orfs

    orfs = find_orfs("CC" + "ATGAAAGGGTTT" + "TAA" + "CC", min_aa=4, include_reverse=False)
    assert len(orfs) == 1 and orfs[0].protein == "MKGF"


# -- 3. structure validation ------------------------------------------------
print("\nStructure validation")


@check("identical sequences align at 100% identity")
def _align():
    from bioseqinsight.core.alignment import compare

    result = compare(UBIQUITIN, UBIQUITIN)
    assert result.identity_percent == 100.0
    assert result.coverage_percent == 100.0


@check("a perfect fragment keeps identity and loses coverage")
def _fragment():
    from bioseqinsight.core.alignment import semi_global

    query = "ACDEFGHIKLMNPQRSTVWY" * 5
    result = semi_global(query, query[20:60])
    assert result.identity_percent == 100.0
    assert abs(result.coverage_percent - 40.0) < 1.0


@check("an unrelated structure is never classified as a match")
def _classify():
    from bioseqinsight.config.settings import Settings
    from bioseqinsight.models.results import MappingLevel
    from bioseqinsight.structures.validation import classify

    settings = Settings()
    assert classify(100.0, 100.0, True, settings) is MappingLevel.M4_EXACT
    # A perfect sequence match without accession confirmation is NOT exact.
    assert classify(100.0, 100.0, False, settings) is MappingLevel.M3_HIGH_CONFIDENCE
    assert classify(5.0, 10.0, False, settings) is MappingLevel.M1_UNVERIFIED


# -- 4. fault tolerance -----------------------------------------------------
print("\nFault tolerance")


@check("a transient failure is retried and recovers")
def _retry():
    from bioseqinsight.services.http_client import HttpClient, RetryPolicy
    from bioseqinsight.services.transport import FakeTransport, ScriptedResponse

    transport = FakeTransport(
        {"example.org": [ScriptedResponse(504, "timeout"), ScriptedResponse(200, "ATOM ok")]}
    )
    client = HttpClient(
        transport=transport,
        policy=RetryPolicy(max_attempts=3, backoff_base_s=0.0, jitter=False),
        sleeper=lambda _s: None,
    )
    outcome = client.get("https://example.org/x")
    assert outcome.ok and outcome.recovered and len(outcome.attempts) == 2


@check("a permanent failure is not retried")
def _no_retry():
    from bioseqinsight.services.http_client import HttpClient, RetryPolicy
    from bioseqinsight.services.transport import FakeTransport, ScriptedResponse

    transport = FakeTransport({"example.org": [ScriptedResponse(404, "missing")]})
    client = HttpClient(
        transport=transport,
        policy=RetryPolicy(max_attempts=3, backoff_base_s=0.0, jitter=False),
        sleeper=lambda _s: None,
    )
    assert len(client.get("https://example.org/x").attempts) == 1


@check("an HTML error page delivered as HTTP 200 is rejected")
def _payload():
    from bioseqinsight.services.http_client import HttpClient, RetryPolicy
    from bioseqinsight.services.transport import FakeTransport, ScriptedResponse

    transport = FakeTransport(
        {
            "example.org": [
                ScriptedResponse(200, "<html>503 backend unavailable</html>"),
                ScriptedResponse(200, "ATOM      1  N   MET A   1"),
            ]
        }
    )
    client = HttpClient(
        transport=transport,
        policy=RetryPolicy(max_attempts=3, backoff_base_s=0.0, jitter=False),
        sleeper=lambda _s: None,
    )
    assert client.get("https://example.org/x", expect_text="ATOM").recovered


# -- 5. end to end ----------------------------------------------------------
print("\nEnd to end")


@check("batch analysis and CSV export")
def _batch():
    from bioseqinsight.config.settings import Settings
    from bioseqinsight.core.alphabet import FastaRecord
    from bioseqinsight.io.export import write_csv
    from bioseqinsight.workflows.batch import BatchRunner

    with tempfile.TemporaryDirectory() as tmp:
        settings = Settings.from_dict(
            {"download_dir": f"{tmp}/s", "cache_dir": f"{tmp}/c",
             "log_file": f"{tmp}/l/a.log", "offline": True}
        )
        rows, summary = BatchRunner(settings).run(
            [FastaRecord("a", "", "ATGAAACCCGGGTTTTAA"), FastaRecord("b", "", UBIQUITIN)]
        )
        assert summary.total == 2 and summary.failed == 0
        path = write_csv(rows, f"{tmp}/out.csv")
        assert Path(path).is_file() and Path(path).stat().st_size > 0


@check("a project round-trips through a zip archive")
def _project():
    from bioseqinsight.config.settings import Settings
    from bioseqinsight.core.alphabet import FastaRecord
    from bioseqinsight.workflows.project import Project

    with tempfile.TemporaryDirectory() as tmp:
        settings = Settings.from_dict(
            {"download_dir": f"{tmp}/s", "cache_dir": f"{tmp}/c",
             "projects_dir": f"{tmp}/p", "log_file": f"{tmp}/l/a.log"}
        )
        project = Project.create(f"{tmp}/proj", name="selfcheck", settings=settings)
        project.add_sequences([FastaRecord("x", "", UBIQUITIN)])
        archive = project.export_zip(f"{tmp}/proj.zip")
        restored = Project.import_zip(archive, f"{tmp}/restored")
        assert restored.manifest.name == "selfcheck"
        assert len(restored.read_sequences()) == 1


@check("results serialise to JSON")
def _json():
    from bioseqinsight.core.sequence import analyze_sequence

    payload = json.loads(json.dumps(analyze_sequence("ATGGCTAGCTAA").to_dict()))
    assert payload["length"] == 12 and payload["software_version"]


@check("the CLI runs")
def _cli():
    import contextlib
    import io

    from bioseqinsight.cli import main

    with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
        code = main(
            ["--quiet", "--offline", "--download-dir", tmp,
             "dna", "ATGAAACCCGGGTTTTAA", "--json", "--out", f"{tmp}/r.json"]
        )
        assert code == 0
        assert json.loads(Path(f"{tmp}/r.json").read_text())["length"] == 18


# -- report -----------------------------------------------------------------
print("\n" + "=" * 60)
total = len(PASSED) + len(FAILED)
if FAILED:
    print(f"{len(PASSED)}/{total} checks passed, {len(FAILED)} FAILED\n")
    for name, error in FAILED:
        print(f"  {name}\n    {error}")
    print("\nRe-run with -v for tracebacks. Please report this with the output of")
    print("'bioseqinsight version' at https://github.com/niloy-999/BioSeqInsight/issues")
    sys.exit(1)

import bioseqinsight  # noqa: E402

print(f"All {total} checks passed. BioSeqInsight {bioseqinsight.__version__} is working.")
sys.exit(0)
