"""Genetic code tables, translation and six-frame translation."""

from __future__ import annotations

import pytest

from bioseqinsight.core.codon import CODON_ORDER, available_tables, get_genetic_code
from bioseqinsight.core.translation import (
    six_frame_translation,
    transcribe,
    translate,
    translate_peptides,
)


class TestGeneticCodeTables:
    def test_every_bundled_table_defines_64_codons(self):
        for table_id in available_tables():
            assert len(get_genetic_code(table_id).forward) == 64

    def test_codon_order_is_the_ncbi_order(self):
        assert CODON_ORDER[0] == "TTT"
        assert CODON_ORDER[-1] == "GGG"
        assert len(set(CODON_ORDER)) == 64

    def test_standard_table_stop_codons(self):
        assert get_genetic_code(1).stops == frozenset({"TAA", "TAG", "TGA"})

    def test_standard_table_start_codons(self):
        assert get_genetic_code(1).starts == frozenset({"ATG", "TTG", "CTG"})

    def test_vertebrate_mitochondrial_differs_where_it_should(self):
        standard, mito = get_genetic_code(1), get_genetic_code(2)
        assert standard.forward["TGA"] == "*" and mito.forward["TGA"] == "W"
        assert standard.forward["AGA"] == "R" and mito.forward["AGA"] == "*"
        assert standard.forward["ATA"] == "I" and mito.forward["ATA"] == "M"

    def test_bacterial_table_has_extra_starts(self):
        assert "GTG" in get_genetic_code(11).starts
        assert "GTG" not in get_genetic_code(1).starts

    def test_unknown_table_is_reported_clearly(self):
        with pytest.raises(KeyError, match="not bundled"):
            get_genetic_code(999)

    def test_tables_are_cached_not_rebuilt(self):
        assert get_genetic_code(1) is get_genetic_code(1)


class TestTranslation:
    @pytest.mark.parametrize(
        "codons,peptide",
        [
            ("ATGGCTTAA", "MA*"),
            ("ATG", "M"),
            ("TGGTGG", "WW"),
            ("AAAAAA", "KK"),
        ],
    )
    def test_known_translations(self, codons, peptide):
        assert translate(codons) == peptide

    def test_incomplete_trailing_codon_is_ignored(self):
        assert translate("ATGGCTA") == "MA"

    def test_to_stop_truncates(self):
        assert translate("ATGGCTTAAGGG", to_stop=True) == "MA"

    def test_ambiguity_codes_translate_to_x(self):
        assert translate("ATGNNN") == "MX"

    def test_frames_shift_the_reading(self):
        sequence = "AATGGCTTAA"
        assert translate(sequence, frame=1) == "MA*"

    def test_invalid_frame_is_rejected(self):
        with pytest.raises(ValueError, match="frame"):
            translate("ATGGCT", frame=3)

    def test_peptides_are_split_on_stops(self):
        assert translate_peptides("ATGGCTTAAATGAAATAA") == ["MA", "MK"]

    def test_empty_input_gives_no_peptides(self):
        assert translate_peptides("") == []

    def test_six_frames_returns_all_six(self):
        frames = six_frame_translation("ATGGCTTAACCCGGGTTTAAA")
        assert set(frames) == {"+1", "+2", "+3", "-1", "-2", "-3"}
        assert frames["+1"].startswith("MA*")

    def test_transcription_matches_translation_input(self):
        assert transcribe("ATGGCT") == "AUGGCU"

    def test_table_choice_changes_the_product(self):
        assert translate("TGA", table=1) == "*"
        assert translate("TGA", table=2) == "W"
