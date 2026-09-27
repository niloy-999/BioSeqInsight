"""Pairwise alignment used for structure-identity validation."""

from __future__ import annotations

import pytest

from bioseqinsight.core.alignment import (
    AlignmentTooLarge,
    compare,
    needleman_wunsch,
    quick_identity,
    semi_global,
)


class TestIdenticalSequences:
    def test_identity_is_one_hundred(self):
        result = compare("ACDEFGHIKL", "ACDEFGHIKL")
        assert result.identity_percent == pytest.approx(100.0)
        assert result.coverage_percent == pytest.approx(100.0)

    def test_fast_path_is_used_for_identical_input(self):
        assert compare("ACDEFGHIKL", "ACDEFGHIKL").mode == "identical"


class TestMismatches:
    def test_single_substitution(self):
        result = needleman_wunsch("ACDEFGHIKL", "ACDEFGHIKM")
        assert result.identity_percent == pytest.approx(90.0)

    def test_completely_different_sequences_score_low(self):
        result = needleman_wunsch("AAAAAAAAAA", "WWWWWWWWWW")
        assert result.identity_percent == pytest.approx(0.0)


class TestFragments:
    def test_a_perfect_fragment_keeps_full_identity(self):
        # This is the case a strict global alignment gets wrong: a PDB chain
        # that is a fragment of the full-length UniProt sequence.
        query = "ACDEFGHIKLMNPQRSTVWY" * 5  # 100 residues
        subject = query[20:60]  # a perfectly matching 40-residue fragment
        result = semi_global(query, subject)
        assert result.identity_percent == pytest.approx(100.0)

    def test_fragment_coverage_reflects_the_missing_part(self):
        query = "ACDEFGHIKLMNPQRSTVWY" * 5
        subject = query[20:60]
        assert semi_global(query, subject).coverage_percent == pytest.approx(40.0)

    def test_global_alignment_penalises_the_same_fragment(self):
        query = "ACDEFGHIKLMNPQRSTVWY" * 5
        subject = query[20:60]
        assert needleman_wunsch(query, subject).identity_percent < 50.0


class TestEdgeCases:
    def test_empty_query(self):
        result = compare("", "ACDEF")
        assert result.identity_percent == 0.0
        assert result.coverage_percent == 0.0

    def test_empty_subject(self):
        assert compare("ACDEF", "").identity_percent == 0.0

    def test_both_empty(self):
        assert compare("", "").identity_percent == 0.0

    def test_case_is_ignored(self):
        assert compare("acdef", "ACDEF").identity_percent == pytest.approx(100.0)

    def test_oversized_alignment_is_refused_rather_than_hanging(self):
        with pytest.raises(AlignmentTooLarge, match="safety limit"):
            needleman_wunsch("A" * 3000, "A" * 3000)


class TestQuickIdentity:
    def test_equal_length_comparison(self):
        assert quick_identity("ACDEF", "ACDEG") == pytest.approx(80.0)

    def test_unequal_lengths_return_zero(self):
        assert quick_identity("ACDEF", "ACDE") == 0.0


class TestInsertions:
    def test_insertion_is_handled(self):
        result = semi_global("ACDEFGHIKL", "ACDEFWGHIKL")
        assert result.identity_percent > 85.0

    def test_deletion_is_handled(self):
        result = semi_global("ACDEFGHIKL", "ACDEGHIKL")
        assert result.identity_percent > 85.0
