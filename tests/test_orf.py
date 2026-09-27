"""Open reading frame detection."""

from __future__ import annotations

import pytest

from bioseqinsight.core.orf import find_orfs, longest_orf, to_forward_coordinates

# ATG AAA GGG TTT TAA  -> M K G F, stop. Padded so the ORF is not at index 0.
SIMPLE = "CC" + "ATGAAAGGGTTT" + "TAA" + "CC"


class TestBasicDetection:
    def test_finds_the_obvious_orf(self):
        orfs = find_orfs(SIMPLE, min_aa=4, include_reverse=False)
        assert len(orfs) == 1
        orf = orfs[0]
        assert orf.protein == "MKGF"
        assert orf.aa_length == 4
        assert orf.start == 2
        assert orf.end == 17
        assert orf.strand == "+"
        assert orf.frame == 3

    def test_min_length_filters(self):
        assert find_orfs(SIMPLE, min_aa=5, include_reverse=False) == []

    def test_no_start_codon_gives_nothing(self):
        assert find_orfs("CCCGGGTTTAAACCC", min_aa=1, include_reverse=False) == []

    def test_results_are_sorted_by_decreasing_length(self):
        sequence = "ATG" + "AAA" * 20 + "TAA" + "ATG" + "AAA" * 5 + "TAA"
        lengths = [orf.aa_length for orf in find_orfs(sequence, min_aa=3, include_reverse=False)]
        assert lengths == sorted(lengths, reverse=True)

    def test_min_aa_must_be_positive(self):
        with pytest.raises(ValueError, match="min_aa"):
            find_orfs(SIMPLE, min_aa=0)


class TestNesting:
    def test_internal_start_does_not_create_a_duplicate_orf(self):
        # Two ATGs, one stop: only the longer ORF should be reported, not a
        # nested copy starting at the second ATG.
        sequence = "ATGAAAATGGGGTAA"
        orfs = find_orfs(sequence, min_aa=2, include_reverse=False)
        forward_frame_one = [o for o in orfs if o.strand == "+" and o.frame == 1]
        assert len(forward_frame_one) == 1
        assert forward_frame_one[0].protein == "MKMG"


class TestStrands:
    def test_reverse_strand_orfs_are_found(self):
        from bioseqinsight.core.codon import reverse_complement

        sequence = reverse_complement("ATG" + "AAA" * 10 + "TAA")
        orfs = find_orfs(sequence, min_aa=5)
        assert any(orf.strand == "-" for orf in orfs)

    def test_reverse_can_be_disabled(self):
        from bioseqinsight.core.codon import reverse_complement

        sequence = reverse_complement("ATG" + "AAA" * 10 + "TAA")
        assert find_orfs(sequence, min_aa=5, include_reverse=False) == []

    def test_forward_coordinate_mapping_is_reversible(self):
        from bioseqinsight.core.codon import reverse_complement

        sequence = reverse_complement("ATG" + "AAA" * 10 + "TAA")
        orf = [o for o in find_orfs(sequence, min_aa=5) if o.strand == "-"][0]
        start, end = to_forward_coordinates(orf, len(sequence))
        assert 0 <= start < end <= len(sequence)
        assert end - start == orf.end - orf.start

    def test_forward_coordinates_are_unchanged_for_plus_strand(self):
        orf = find_orfs(SIMPLE, min_aa=4, include_reverse=False)[0]
        assert to_forward_coordinates(orf, len(SIMPLE)) == (orf.start, orf.end)


class TestUnterminated:
    def test_orf_running_off_the_end_is_flagged(self):
        sequence = "ATG" + "AAA" * 12  # no stop codon
        orfs = [o for o in find_orfs(sequence, min_aa=5, include_reverse=False)]
        assert orfs and orfs[0].unterminated is True

    def test_unterminated_orfs_can_be_excluded(self):
        sequence = "ATG" + "AAA" * 12
        assert find_orfs(
            sequence, min_aa=5, include_reverse=False, include_unterminated=False
        ) == []


class TestAlternativeStarts:
    def test_atg_only_by_default(self):
        sequence = "GTG" + "AAA" * 10 + "TAA"
        assert find_orfs(sequence, min_aa=5, include_reverse=False, table=11) == []

    def test_alternative_starts_can_be_enabled(self):
        sequence = "GTG" + "AAA" * 10 + "TAA"
        orfs = find_orfs(
            sequence, min_aa=5, include_reverse=False, table=11, require_atg=False
        )
        assert orfs and orfs[0].protein.startswith("M")

    def test_longest_orf_helper(self):
        assert longest_orf(SIMPLE, min_aa=4, include_reverse=False).protein == "MKGF"

    def test_longest_orf_returns_none_when_absent(self):
        assert longest_orf("CCCCCC", min_aa=5) is None
