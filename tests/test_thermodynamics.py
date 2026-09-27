"""Melting-temperature calculation."""

from __future__ import annotations

import pytest

from bioseqinsight.core.thermodynamics import (
    NN_PARAMS,
    ThermodynamicsError,
    is_self_complementary,
    melting_temperature,
    nearest_neighbor_thermodynamics,
    nearest_neighbor_tm,
    salt_correction,
    wallace_tm,
)


class TestWallace:
    @pytest.mark.parametrize(
        "sequence,expected",
        [("AAAA", 8.0), ("GGGG", 16.0), ("ATGC", 12.0), ("ATGCATGCATGCATGCATGC", 60.0)],
    )
    def test_known_values(self, sequence, expected):
        assert wallace_tm(sequence) == pytest.approx(expected)

    def test_empty_sequence_is_rejected(self):
        with pytest.raises(ThermodynamicsError):
            wallace_tm("")


class TestNearestNeighbour:
    def test_parameter_table_is_symmetric(self):
        # Each doublet and its reverse complement must share parameters.
        complement = str.maketrans("ACGT", "TGCA")
        for doublet, values in NN_PARAMS.items():
            mirrored = doublet.translate(complement)[::-1]
            assert NN_PARAMS[mirrored] == values

    def test_enthalpy_is_negative_for_any_duplex(self):
        delta_h, delta_s = nearest_neighbor_thermodynamics("ATGCATGCATGC")
        assert delta_h < 0
        assert delta_s < 0

    def test_too_short_for_nearest_neighbour(self):
        with pytest.raises(ThermodynamicsError, match="at least 2"):
            nearest_neighbor_thermodynamics("A")

    def test_gc_rich_melts_higher_than_at_rich(self):
        gc_rich = nearest_neighbor_tm("GCGCGCGCGCGCGCGCGCGC")
        at_rich = nearest_neighbor_tm("ATATATATATATATATATAT")
        assert gc_rich > at_rich + 20

    def test_values_are_physically_plausible(self):
        # A 20-mer at 50% GC should melt somewhere in the 50-65 C range.
        tm = nearest_neighbor_tm("ATGCATGCATGCATGCATGC")
        assert 45.0 < tm < 70.0

    def test_longer_duplex_melts_higher(self):
        short = nearest_neighbor_tm("ATGCATGCATGCATGC")
        long = nearest_neighbor_tm("ATGCATGCATGCATGCATGCATGCATGC")
        assert long > short

    def test_higher_salt_raises_tm(self):
        low = nearest_neighbor_tm("ATGCATGCATGCATGCATGC", na_mM=10)
        high = nearest_neighbor_tm("ATGCATGCATGCATGCATGC", na_mM=1000)
        assert high > low

    def test_zero_salt_is_rejected(self):
        with pytest.raises(ThermodynamicsError, match="Sodium"):
            salt_correction(-100.0, 20, 0.0)

    def test_zero_primer_concentration_is_rejected(self):
        with pytest.raises(ThermodynamicsError, match="Primer"):
            nearest_neighbor_tm("ATGCATGCATGC", primer_nM=0)

    def test_self_complementary_detection(self):
        assert is_self_complementary("GAATTC")
        assert not is_self_complementary("GAATTA")


class TestMethodSelection:
    def test_short_oligo_uses_wallace(self):
        result = melting_temperature("ATGCATGCATGCATGC")
        assert result.method == "wallace"

    def test_long_sequence_uses_nearest_neighbour(self):
        result = melting_temperature("ATGCATGCATGCATGCATGCATGC")
        assert result.method == "nearest_neighbor"

    def test_wallace_on_a_gene_carries_a_warning(self):
        result = melting_temperature("ATGC" * 100, method="wallace")
        assert result.note and "14-20" in result.note

    def test_very_long_nearest_neighbour_carries_a_warning(self):
        result = melting_temperature("ATGC" * 50)
        assert result.note and "two-state" in result.note.lower()

    def test_wallace_value_is_always_reported_alongside(self):
        result = melting_temperature("ATGCATGCATGCATGCATGCATGC")
        assert result.wallace_c == pytest.approx(72.0)

    def test_unknown_method_is_rejected(self):
        with pytest.raises(ThermodynamicsError, match="Unknown"):
            melting_temperature("ATGCATGC", method="magic")

    def test_no_canonical_bases_is_rejected(self):
        with pytest.raises(ThermodynamicsError, match="undefined"):
            melting_temperature("NNNNNNNN")

    def test_result_records_the_conditions(self):
        result = melting_temperature("ATGCATGCATGCATGCATGCATGC", primer_nM=500, na_mM=100)
        assert result.primer_nM == 500
        assert result.na_mM == 100
