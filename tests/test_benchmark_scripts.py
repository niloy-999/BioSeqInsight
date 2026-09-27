"""Sanity checks on the benchmark scripts themselves.

The scripts under benchmarks/scripts/ are not part of the installed package,
but they produce the numbers cited in the manuscript, so their arithmetic
gets the same scrutiny as anything in src/. This module was added after a
live run showed a per-provider table where retrieved + not_found + timeouts
+ service_errors did not sum to requests -- not because any number was
wrong, but because ESM Atlas's length-limit rejections (status
"invalid_input") were not counted anywhere in the summary at all.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from run_structure_benchmark import summarise  # noqa: E402


def _row(accession, provider, status, mapping="M0", recovered=False, response_time_s=None):
    return {
        "accession": accession,
        "provider": provider,
        "status": status,
        "mapping_level": mapping,
        "recovered": recovered,
        "response_time_s": response_time_s,
    }


class TestSummariseAccounting:
    """Every request must land in exactly one reported bucket."""

    def test_every_status_category_is_counted(self):
        rows = [
            _row("A", "esmatlas", "ok", mapping="M4", response_time_s=1.0),
            _row("B", "esmatlas", "not_found"),
            _row("C", "esmatlas", "timeout"),
            _row("D", "esmatlas", "service_error"),
            _row("E", "esmatlas", "invalid_input"),
            _row("F", "esmatlas", "skipped"),
        ]
        resolved = [{"accession": a, "stratum": "short"} for a in "ABCDEF"]
        summary = summarise(rows, resolved, [])["per_provider"]["esmatlas"]
        assert summary["requests"] == 6
        assert summary["accounting_check"] == "ok"
        assert (
            summary["retrieved"] + summary["not_found"] + summary["timeouts"]
            + summary["service_errors"] + summary["invalid_input"] + summary["skipped"]
        ) == summary["requests"]

    def test_invalid_input_is_reported_not_silently_dropped(self):
        # This is the exact case a live run hit: sequences over ESM Atlas's
        # length limit are refused before any HTTP call and must still
        # appear somewhere in the summary, not just vanish from the totals.
        rows = [_row("A", "esmatlas", "invalid_input") for _ in range(23)] + [
            _row("B", "esmatlas", "ok", mapping="M4", response_time_s=0.5)
        ]
        resolved = [{"accession": "A", "stratum": "long"}, {"accession": "B", "stratum": "short"}]
        summary = summarise(rows, resolved, [])["per_provider"]["esmatlas"]
        assert summary["invalid_input"] == 23
        assert summary["requests"] == 24
        assert summary["accounting_check"] == "ok"

    def test_accounting_check_flags_a_genuine_mismatch(self):
        # A status value the summariser does not know about must be visible
        # as a mismatch, not silently ignored.
        rows = [_row("A", "esmatlas", "some_new_status_nobody_added_yet")]
        resolved = [{"accession": "A", "stratum": "short"}]
        summary = summarise(rows, resolved, [])["per_provider"]["esmatlas"]
        assert summary["accounting_check"] != "ok"
        assert "MISMATCH" in summary["accounting_check"]

    def test_median_response_time_ignores_non_numeric_entries(self):
        rows = [
            _row("A", "rcsb", "ok", mapping="M4", response_time_s=2.0),
            _row("B", "rcsb", "ok", mapping="M4", response_time_s=4.0),
            _row("C", "rcsb", "not_found", response_time_s=None),
        ]
        resolved = [{"accession": a, "stratum": "short"} for a in "ABC"]
        summary = summarise(rows, resolved, [])["per_provider"]["rcsb"]
        assert summary["median_response_time_s"] == 3.0

    def test_empty_provider_subset_does_not_divide_by_zero(self):
        summary = summarise([], [], [])["per_provider"]
        assert summary == {}


class TestBestMappingDistribution:
    """Regression tests for a real defect: a protein that got no structure
    from any provider vanished from best_mapping_distribution instead of
    being counted under M0, so the distribution silently summed to one less
    than the number of resolved candidates. Caught by a live run where
    M0+M1+M2+M3+M4 summed to 68 against 69 resolved proteins.
    """

    def test_a_protein_with_no_structure_anywhere_is_counted_as_m0(self):
        rows = [
            _row("A", "rcsb", "not_found"),
            _row("A", "alphafold", "not_found"),
            _row("A", "esmatlas", "service_error"),
            _row("B", "rcsb", "ok", mapping="M4"),
        ]
        resolved = [{"accession": "A", "stratum": "short"}, {"accession": "B", "stratum": "short"}]
        summary = summarise(rows, resolved, [])
        assert summary["best_mapping_distribution"]["M0"] == 1
        assert summary["best_mapping_distribution"]["M4"] == 1

    def test_distribution_always_sums_to_the_resolved_count(self):
        rows = [
            _row("A", "rcsb", "not_found"),
            _row("B", "rcsb", "ok", mapping="M2"),
            _row("C", "rcsb", "ok", mapping="M4"),
            _row("D", "rcsb", "service_error"),
        ]
        resolved = [{"accession": a, "stratum": "short"} for a in "ABCD"]
        summary = summarise(rows, resolved, [])
        assert sum(summary["best_mapping_distribution"].values()) == len(resolved)

    def test_a_resolved_protein_absent_from_rows_entirely_is_still_m0(self):
        # Defensive: even if a provider never produced a row for some
        # resolved accession at all (rather than an explicit not_found row),
        # it must still be counted, not silently dropped.
        rows = [_row("A", "rcsb", "ok", mapping="M4")]
        resolved = [{"accession": "A", "stratum": "short"}, {"accession": "B", "stratum": "short"}]
        summary = summarise(rows, resolved, [])
        assert summary["best_mapping_distribution"]["M0"] == 1
        assert sum(summary["best_mapping_distribution"].values()) == 2

    def test_rate_percentages_are_unaffected_by_the_fix(self):
        # any_structure_rate_percent and exact_mapping_rate_percent are
        # derived without relying on distribution["M0"], so they were
        # already correct even while the distribution table was not; this
        # pins that they remain correct after the fix too.
        rows = [
            _row("A", "rcsb", "not_found"),
            _row("B", "rcsb", "ok", mapping="M3"),
            _row("C", "rcsb", "ok", mapping="M4"),
        ]
        resolved = [{"accession": a, "stratum": "short"} for a in "ABC"]
        summary = summarise(rows, resolved, [])
        assert summary["exact_mapping_rate_percent"] == pytest.approx(100.0 / 3, abs=0.1)
        assert summary["any_structure_rate_percent"] == pytest.approx(200.0 / 3, abs=0.1)


class TestAttemptedSuccessRate:
    """attempted_success_rate_percent excludes client-side refusals from the
    denominator, separating "this provider has no coverage / we didn't even
    ask" from "we asked and it failed" -- the distinction that made ESM
    Atlas's raw 66.7% figure easy to misread as unreliability when it was
    actually the length-limit exclusion rate.
    """

    def test_excludes_invalid_input_from_the_denominator(self):
        rows = [_row(f"L{i}", "esmatlas", "invalid_input") for i in range(23)] + [
            _row(f"S{i}", "esmatlas", "ok", mapping="M4", response_time_s=0.5) for i in range(46)
        ]
        resolved = [{"accession": r["accession"], "stratum": "short"} for r in rows]
        summary = summarise(rows, resolved, [])["per_provider"]["esmatlas"]
        assert summary["success_rate_percent"] == pytest.approx(100.0 * 46 / 69, abs=0.1)
        assert summary["attempted_requests"] == 46
        assert summary["attempted_success_rate_percent"] == 100.0

    def test_none_when_nothing_was_attempted(self):
        rows = [_row("A", "esmatlas", "invalid_input")]
        resolved = [{"accession": "A", "stratum": "long"}]
        summary = summarise(rows, resolved, [])["per_provider"]["esmatlas"]
        assert summary["attempted_requests"] == 0
        assert summary["attempted_success_rate_percent"] is None
