"""Nucleotide analysis: composition, GC, transcription, reverse complement."""

from __future__ import annotations

import pytest

from bioseqinsight.core.alphabet import SequenceError
from bioseqinsight.core.codon import complement, reverse_complement
from bioseqinsight.core.sequence import (
    SequenceAnalyzer,
    analyze_sequence,
    gc_content,
    nucleotide_composition,
)


class TestGCContent:
    @pytest.mark.parametrize(
        "sequence,expected",
        [("GGCC", 100.0), ("AATT", 0.0), ("ATGC", 50.0), ("GGGGAATT", 50.0)],
    )
    def test_known_values(self, sequence, expected):
        assert gc_content(sequence) == pytest.approx(expected)

    def test_ambiguity_codes_are_excluded_from_both_terms(self):
        # NNNN must not depress GC%: 2 of 4 real bases are G/C.
        assert gc_content("ATGCNNNN") == pytest.approx(50.0)

    def test_empty_sequence_is_zero_not_an_error(self):
        assert gc_content("") == 0.0

    def test_case_insensitive(self):
        assert gc_content("gcgc") == pytest.approx(100.0)


class TestComposition:
    def test_counts_are_exact(self):
        composition = nucleotide_composition("AATTTGGGGCC")
        assert composition == {"A": 2, "T": 3, "G": 4, "C": 2, "N": 0}

    def test_ambiguity_codes_collapse_into_n(self):
        composition = nucleotide_composition("ATGCNRYW")
        assert composition["N"] == 4

    def test_percentages_sum_to_one_hundred(self):
        result = analyze_sequence("ATGCATGCATGC")
        assert sum(result.composition_percent.values()) == pytest.approx(100.0)


class TestComplement:
    def test_reverse_complement_of_known_sequence(self):
        assert reverse_complement("ATGC") == "GCAT"

    def test_double_reverse_complement_is_identity(self):
        sequence = "ATGCGATCGTAGCTAGCTA"
        assert reverse_complement(reverse_complement(sequence)) == sequence

    def test_ambiguity_codes_are_complemented_correctly(self):
        assert complement("RYSWKM") == "YRSWMK"
        assert complement("N") == "N"

    def test_palindrome_is_its_own_reverse_complement(self):
        assert reverse_complement("GAATTC") == "GAATTC"


class TestTranscription:
    def test_thymine_becomes_uracil(self):
        assert SequenceAnalyzer("ATGCTT").transcribe() == "AUGCUU"

    def test_rna_input_is_accepted_and_normalised(self):
        assert SequenceAnalyzer("AUGCUU").sequence == "ATGCTT"


class TestAnalyzer:
    def test_length_excludes_stripped_characters(self):
        analyzer = SequenceAnalyzer(">header\nATG CTA\n123")
        assert analyzer.get_length() == 6

    def test_empty_input_raises(self):
        with pytest.raises(SequenceError):
            SequenceAnalyzer("")

    def test_invalid_input_raises(self):
        with pytest.raises(SequenceError):
            SequenceAnalyzer("ATGC@@")

    def test_analyze_returns_a_populated_result(self):
        result = analyze_sequence("ATG" + "GCT" * 20 + "TAA", identifier="demo")
        assert result.identifier == "demo"
        assert result.length == 66
        assert result.software_version
        assert result.tm is not None
        assert result.reverse_complement is not None

    def test_include_sequences_false_omits_bulky_strings(self):
        result = SequenceAnalyzer("ATG" + "GCT" * 20 + "TAA").analyze(include_sequences=False)
        assert result.reverse_complement is None
        assert result.mrna is None
        assert result.peptides == []

    def test_ambiguity_codes_produce_a_warning(self):
        result = analyze_sequence("ATGCNNNNATGC")
        assert result.ambiguous_count == 4
        assert any("ambiguity" in warning for warning in result.warnings)

    def test_result_is_json_serialisable(self):
        import json

        result = analyze_sequence("ATGGCTAGCTAA")
        assert json.loads(json.dumps(result.to_dict()))["length"] == 12

    def test_gc_skew_windows(self):
        skew = SequenceAnalyzer("GGGGCCCC" * 4).get_gc_skew(window=8)
        assert len(skew) == 4
        assert all(value == pytest.approx(0.0) for _, value in skew)
