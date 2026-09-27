"""Protein physicochemical descriptors."""

from __future__ import annotations

import pytest

from bioseqinsight.core.alphabet import SequenceError
from bioseqinsight.core.protein import (
    AVERAGE_RESIDUE_MASS,
    KYTE_DOOLITTLE,
    ProteinAnalyzer,
    aliphatic_index,
    amino_acid_composition,
    aromaticity,
    charge_at_ph,
    extinction_coefficient,
    gravy,
    hydropathy_profile,
    isoelectric_point,
    molecular_weight,
    secondary_structure_sketch,
)

UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)


class TestMolecularWeight:
    def test_single_residue_equals_its_average_mass(self):
        for residue, mass in AVERAGE_RESIDUE_MASS.items():
            assert molecular_weight(residue) == pytest.approx(mass, abs=1e-6)

    def test_ubiquitin_matches_the_published_value(self):
        # Ubiquitin is 8.6 kDa; the accepted average mass is 8564.8 Da.
        assert molecular_weight(UBIQUITIN) == pytest.approx(8564.8, abs=0.5)

    def test_water_is_removed_once_per_peptide_bond(self):
        two = molecular_weight("AA")
        assert two == pytest.approx(2 * AVERAGE_RESIDUE_MASS["A"] - 18.0153, abs=1e-6)

    def test_ambiguity_codes_are_excluded(self):
        assert molecular_weight("AAX") == pytest.approx(molecular_weight("AA"), abs=1e-6)

    def test_only_ambiguity_codes_is_an_error(self):
        with pytest.raises(SequenceError, match="no standard amino acids"):
            molecular_weight("XXXX")


class TestHydropathy:
    def test_gravy_of_a_single_residue(self):
        for residue, value in KYTE_DOOLITTLE.items():
            assert gravy(residue) == pytest.approx(value)

    def test_poly_isoleucine_is_maximally_hydrophobic(self):
        assert gravy("I" * 20) == pytest.approx(4.5)

    def test_poly_arginine_is_maximally_hydrophilic(self):
        assert gravy("R" * 20) == pytest.approx(-4.5)

    def test_ubiquitin_is_hydrophilic_overall(self):
        assert -1.0 < gravy(UBIQUITIN) < 0.0

    def test_profile_length(self):
        profile = hydropathy_profile("A" * 20, window=9)
        assert len(profile) == 12

    def test_profile_is_empty_when_sequence_is_shorter_than_window(self):
        assert hydropathy_profile("AAA", window=9) == []

    def test_window_must_be_positive(self):
        with pytest.raises(ValueError, match="window"):
            hydropathy_profile("AAAA", window=0)


class TestChargeAndPI:
    def test_charge_falls_monotonically_with_ph(self):
        values = [charge_at_ph(UBIQUITIN, ph) for ph in (3, 5, 7, 9, 11)]
        assert values == sorted(values, reverse=True)

    def test_pi_is_where_charge_is_zero(self):
        pi = isoelectric_point(UBIQUITIN)
        assert abs(charge_at_ph(UBIQUITIN, pi)) < 0.05

    def test_acidic_protein_has_low_pi(self):
        assert isoelectric_point("DDDDEEEEDDDD") < 4.5

    def test_basic_protein_has_high_pi(self):
        assert isoelectric_point("KKKKRRRRKKKK") > 10.0

    def test_pi_is_within_the_scale(self):
        assert 0.0 <= isoelectric_point(UBIQUITIN) <= 14.0


class TestOtherDescriptors:
    def test_aromaticity_counts_fwy(self):
        assert aromaticity("FWYAAAAAAA") == pytest.approx(0.3)

    def test_aromaticity_of_a_non_aromatic_peptide_is_zero(self):
        assert aromaticity("AAAAGGGG") == 0.0

    def test_extinction_coefficient_uses_gill_von_hippel(self):
        reduced, cystines = extinction_coefficient("WWYYCC")
        assert reduced == 2 * 5500 + 2 * 1490
        assert cystines == reduced + 125

    def test_extinction_without_aromatics_is_zero(self):
        assert extinction_coefficient("AAAA") == (0, 0)

    def test_aliphatic_index_of_poly_alanine(self):
        assert aliphatic_index("A" * 10) == pytest.approx(100.0)

    def test_composition_covers_all_twenty_residues(self):
        composition = amino_acid_composition("ACDEF")
        assert len(composition) == 20
        assert composition["A"] == 1 and composition["G"] == 0


class TestSecondaryStructureSketch:
    def test_percentages_sum_to_one_hundred(self):
        sketch = secondary_structure_sketch(UBIQUITIN)
        total = sketch.helix_percent + sketch.sheet_percent + sketch.coil_percent
        assert total == pytest.approx(100.0)

    def test_visual_length_matches_the_sequence(self):
        sketch = secondary_structure_sketch(UBIQUITIN)
        assert len(sketch.visual) == len(UBIQUITIN)

    def test_method_string_says_it_is_not_a_prediction(self):
        sketch = secondary_structure_sketch(UBIQUITIN)
        assert "sketch" in sketch.method or "illustrative" in sketch.method

    def test_poly_glutamate_is_called_helical(self):
        # E has the highest helix propensity in the table.
        assert secondary_structure_sketch("E" * 30).helix_percent == pytest.approx(100.0)


class TestAnalyzer:
    def test_full_analysis_is_populated(self):
        result = ProteinAnalyzer(UBIQUITIN, identifier="ubq").analyze()
        assert result.identifier == "ubq"
        assert result.length == 76
        assert result.molecular_weight > 8000
        assert result.secondary_structure is not None

    def test_sketch_can_be_skipped(self):
        assert ProteinAnalyzer(UBIQUITIN).analyze(include_sketch=False).secondary_structure is None

    def test_ambiguity_codes_produce_a_warning(self):
        result = ProteinAnalyzer("MKTAYXXX").analyze()
        assert result.unknown_residues == 3
        assert any("ambiguity" in warning for warning in result.warnings)

    def test_empty_input_is_rejected(self):
        with pytest.raises(SequenceError):
            ProteinAnalyzer("")

    def test_composition_percentages_sum_to_one_hundred_for_standard_residues(self):
        result = ProteinAnalyzer(UBIQUITIN).analyze()
        assert sum(result.composition_percent.values()) == pytest.approx(100.0)

    def test_result_is_json_serialisable(self):
        import json

        result = ProteinAnalyzer(UBIQUITIN).analyze()
        assert json.loads(json.dumps(result.to_dict()))["length"] == 76
