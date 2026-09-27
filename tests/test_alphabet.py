"""Input normalisation, validation and FASTA handling."""

from __future__ import annotations

import pytest

from bioseqinsight.core.alphabet import (
    FastaRecord,
    SequenceError,
    canonical_dna,
    clean_dna,
    clean_protein,
    detect_alphabet,
    invalid_characters,
    iter_fasta_file,
    looks_like_dna,
    looks_like_protein,
    parse_fasta,
    strip_headers_and_whitespace,
    validate_dna,
    validate_protein,
    write_fasta,
)


class TestCleaning:
    def test_whitespace_and_digits_are_removed(self):
        assert clean_dna("  atg cta\n123 gca ") == "ATGCTAGCA"

    def test_lowercase_is_accepted(self):
        assert clean_dna("atgc") == "ATGC"

    def test_uracil_is_folded_to_thymine(self):
        assert clean_dna("AUGC") == "ATGC"

    def test_fasta_header_is_dropped(self):
        assert clean_dna(">seq1 description\nATGC\nGGTA") == "ATGCGGTA"

    def test_comment_lines_are_dropped(self):
        assert clean_dna("; a comment\nATGC") == "ATGC"

    def test_ambiguity_codes_survive_cleaning_but_not_canonicalisation(self):
        assert clean_dna("ATGCNRY") == "ATGCNRY"
        assert canonical_dna("ATGCNRY") == "ATGC"

    def test_empty_input_gives_empty_string(self):
        assert clean_dna("") == ""
        assert clean_protein(None or "") == ""

    def test_protein_ambiguity_codes_are_kept(self):
        assert clean_protein("MKXZB") == "MKXZB"
        assert clean_protein("MKXZB", keep_ambiguous=False) == "MK"


class TestValidation:
    def test_valid_dna_round_trips(self):
        assert validate_dna(" atg\ncta ") == "ATGCTA"

    def test_empty_input_is_rejected_with_guidance(self):
        with pytest.raises(SequenceError, match="Empty input"):
            validate_dna("   ")

    def test_invalid_dna_character_is_named(self):
        with pytest.raises(SequenceError, match="'Z'"):
            validate_dna("ATGZ")

    def test_symbols_are_rejected_but_position_numbers_are_not(self):
        # Digits and gaps appear in sequences pasted out of database records
        # and are ignored; a genuine stray symbol is reported.
        assert validate_protein("1 MKTAY 5") == "MKTAY"
        with pytest.raises(SequenceError, match="'@'"):
            validate_protein("MKT@AY")

    def test_valid_protein_round_trips(self):
        assert validate_protein("mkt ay") == "MKTAY"

    def test_invalid_characters_lists_each_offender_once(self):
        assert invalid_characters("ATGZZQ", "ACGT") == ["Z", "Q"]


class TestDetection:
    def test_dna_is_detected(self):
        assert looks_like_dna("ATGCGATCGATCGATC")
        assert detect_alphabet("ATGCGATCGATCGATC") == "dna"

    def test_protein_is_detected(self):
        assert looks_like_protein("MKTAYIAKQRQISFVKSHFSRQ")
        assert detect_alphabet("MKTAYIAKQRQISFVKSHFSRQ") == "protein"

    def test_a_peptide_of_acgt_residues_is_not_called_dna_when_it_has_others(self):
        # ACGT are all valid amino acids too; the discriminator is the presence
        # of residues that cannot be nucleotides.
        assert detect_alphabet("ACGTACGTMKWL") == "protein"

    def test_very_short_input_is_unknown(self):
        assert detect_alphabet("AT") == "unknown"

    def test_gibberish_is_unknown(self):
        assert detect_alphabet("!!!@@@###") == "unknown"


class TestFasta:
    def test_single_record(self):
        records = parse_fasta(">seq1 my description\nATGC\nGGTA")
        assert len(records) == 1
        assert records[0].identifier == "seq1"
        assert records[0].description == "my description"
        assert records[0].sequence == "ATGCGGTA"

    def test_multiple_records(self):
        records = parse_fasta(">a\nATG\n>b\nCCC\n>c\nGGG")
        assert [r.identifier for r in records] == ["a", "b", "c"]
        assert [r.sequence for r in records] == ["ATG", "CCC", "GGG"]

    def test_headerless_text_becomes_one_record(self):
        records = parse_fasta("ATGCATGC")
        assert len(records) == 1
        assert records[0].identifier == "sequence_1"

    def test_empty_records_are_skipped(self):
        records = parse_fasta(">empty\n\n>real\nATG")
        assert [r.identifier for r in records] == ["real"]

    def test_blank_input_gives_no_records(self):
        assert parse_fasta("") == []

    def test_write_fasta_wraps_lines(self):
        text = write_fasta([FastaRecord("s", "", "A" * 130)], line_width=60)
        lines = text.strip().splitlines()
        assert lines[0] == ">s"
        assert [len(line) for line in lines[1:]] == [60, 60, 10]

    def test_round_trip(self):
        original = [FastaRecord("x", "desc here", "ATGCATGC")]
        assert parse_fasta(write_fasta(original))[0].sequence == "ATGCATGC"

    def test_file_streaming_matches_text_parsing(self, tmp_path):
        path = tmp_path / "multi.fasta"
        path.write_text(">a\nATG\n>b\nCCCGGG\n", encoding="utf-8")
        streamed = list(iter_fasta_file(str(path)))
        assert [(r.identifier, r.sequence) for r in streamed] == [("a", "ATG"), ("b", "CCCGGG")]

    def test_strip_headers_keeps_letters_only(self):
        assert strip_headers_and_whitespace(">h\nAT-G*C 12\n") == "ATGC"
