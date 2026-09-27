"""Motif searching."""

from __future__ import annotations

import pytest

from bioseqinsight.core.motifs import (
    MotifError,
    find_motif,
    gc_skew,
    iupac_to_regex,
    restriction_summary,
)


class TestIUPACExpansion:
    @pytest.mark.parametrize(
        "motif,pattern",
        [("ACGT", "ACGT"), ("N", "[ACGT]"), ("R", "[AG]"), ("WGATAR", "[AT]GATA[AG]")],
    )
    def test_expansion(self, motif, pattern):
        assert iupac_to_regex(motif) == pattern

    def test_empty_motif_is_rejected(self):
        with pytest.raises(MotifError, match="Empty motif"):
            iupac_to_regex("")


class TestSearch:
    def test_exact_match_positions_are_zero_based(self):
        hits = find_motif("AAGAATTCAA", "GAATTC")
        assert [(h.start, h.end) for h in hits] == [(2, 8)]

    def test_absent_motif_gives_no_hits(self):
        assert find_motif("AAAAAAA", "GAATTC") == []

    def test_overlapping_occurrences_are_all_reported(self):
        # AAAA contains three overlapping AA motifs.
        assert len(find_motif("AAAA", "AA")) == 3

    def test_degenerate_motif_matches_every_variant(self):
        # WGATAR = [AT]GATA[AG]: TGATAA and AGATAG both match, CGATAG does not.
        hits = find_motif("TGATAAccAGATAGccCGATAG", "WGATAR")
        assert [(h.start, h.end) for h in hits] == [(0, 6), (8, 14)]

    def test_both_strands_reports_forward_coordinates(self):
        sequence = "AAAAGAATTCAAAA"
        hits = find_motif(sequence, "GAATTC", both_strands=True)
        # EcoRI is palindromic, so it is found once on each strand at the
        # same forward coordinates.
        assert len(hits) == 2
        assert {(h.start, h.end) for h in hits} == {(4, 10)}
        assert {h.strand for h in hits} == {"+", "-"}

    def test_non_palindromic_motif_on_the_reverse_strand(self):
        # "GGGGG" reverse complement is "CCCCC".
        hits = find_motif("AACCCCCAA", "GGGGG", both_strands=True)
        assert [h.strand for h in hits] == ["-"]
        assert hits[0].start == 2

    def test_ambiguity_in_the_subject_does_not_shift_coordinates(self):
        hits = find_motif("NNGAATTC", "GAATTC")
        assert hits[0].start == 2

    def test_regex_mode(self):
        hits = find_motif("ATGAAATTTATG", "ATG", regex=True)
        assert len(hits) == 2

    def test_invalid_regex_is_reported(self):
        with pytest.raises(MotifError, match="Invalid pattern"):
            find_motif("ATGC", "([", regex=True)

    def test_invalid_iupac_code_is_reported(self):
        with pytest.raises(MotifError):
            find_motif("ATGC", "Z")

    def test_empty_sequence_gives_no_hits(self):
        assert find_motif("", "ATG") == []


class TestGCSkew:
    def test_balanced_windows_are_zero(self):
        assert all(value == 0.0 for _, value in gc_skew("GGGCCC" * 4, window=6))

    def test_g_rich_window_is_positive(self):
        assert gc_skew("GGGGGGGGGG", window=10)[0][1] == pytest.approx(1.0)

    def test_window_must_be_positive(self):
        with pytest.raises(ValueError, match="window"):
            gc_skew("ATGC", window=0)


class TestRestrictionSummary:
    def test_user_supplied_sites_are_located(self):
        summary = restriction_summary("AAGAATTCAA", {"EcoRI": "GAATTC"})
        assert summary["EcoRI"] == [2]

    def test_invalid_sites_are_skipped_silently(self):
        assert restriction_summary("ATGC", {"Bad": "ZZZ"}) == {}
