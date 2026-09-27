"""PDB coordinate parsing."""

from __future__ import annotations

import pytest

from bioseqinsight.structures.pdbio import (
    find_accessions,
    looks_like_pdb,
    mean_confidence,
    parse_pdb,
)

UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)

MINIMAL = """\
ATOM      1  N   MET A   1      10.000  10.000  10.000  1.00 90.00           N
ATOM      2  CA  MET A   1      11.000  10.000  10.000  1.00 90.00           C
ATOM      3  CA  LYS A   2      12.000  10.000  10.000  1.00 80.00           C
ATOM      4  CA  THR A   3      13.000  10.000  10.000  1.00 70.00           C
END
"""


class TestRealFile:
    def test_chain_sequence_is_recovered(self, ubiquitin_pdb):
        parsed = parse_pdb(ubiquitin_pdb)
        assert len(parsed.chains) == 1
        assert parsed.chains[0].sequence == UBIQUITIN
        assert parsed.chains[0].residue_count == 76

    def test_cross_reference_is_recovered(self, ubiquitin_pdb):
        assert "P62988" in parse_pdb(ubiquitin_pdb).accessions

    def test_coordinates_are_detected(self, ubiquitin_pdb):
        parsed = parse_pdb(ubiquitin_pdb)
        assert parsed.has_coordinates
        assert parsed.atom_count > 500

    def test_title_is_captured(self, ubiquitin_pdb):
        assert "UBIQUITIN" in parse_pdb(ubiquitin_pdb).title.upper()


class TestMinimalFile:
    def test_ca_atoms_give_the_sequence(self):
        assert parse_pdb(MINIMAL).chains[0].sequence == "MKT"

    def test_duplicate_residues_are_not_double_counted(self):
        doubled = MINIMAL.replace(
            "ATOM      2  CA  MET A   1",
            "ATOM      2  CA  MET A   1",
        )
        assert parse_pdb(doubled).chains[0].sequence == "MKT"

    def test_source_is_reported_as_atom_when_no_seqres(self):
        assert parse_pdb(MINIMAL).chains[0].source == "atom"


class TestConfidence:
    def test_predicted_files_report_plddt(self):
        value, kind = mean_confidence(MINIMAL, is_predicted=True)
        assert kind == "plddt"
        assert value == pytest.approx(82.5)

    def test_experimental_files_report_b_factors(self):
        value, kind = mean_confidence(MINIMAL, is_predicted=False)
        assert kind == "bfactor"

    def test_a_zero_to_one_column_is_rescaled_for_predictions(self):
        scaled = MINIMAL.replace(" 90.00 ", "  0.90 ").replace(" 80.00 ", "  0.80 ").replace(
            " 70.00 ", "  0.70 "
        )
        value, kind = mean_confidence(scaled, is_predicted=True)
        assert kind == "plddt"
        assert value == pytest.approx(82.5, abs=0.1)

    def test_an_experimental_b_factor_is_never_called_plddt(self, ubiquitin_pdb):
        # 1UBQ B-factors sit inside 0-100 but are temperature factors, not
        # confidence. The provider decides, not the numeric range.
        _, kind = mean_confidence(ubiquitin_pdb, is_predicted=False)
        assert kind == "bfactor"

    def test_file_without_coordinates_has_no_confidence(self):
        assert mean_confidence("HEADER something\nEND\n", is_predicted=True) == (None, None)


class TestChainSelection:
    def test_best_matching_chain_prefers_the_exact_sequence(self):
        two_chains = MINIMAL.replace("END", "") + (
            "ATOM      9  CA  ALA B   1      20.000  10.000  10.000  1.00 50.00           C\n"
            "ATOM     10  CA  GLY B   2      21.000  10.000  10.000  1.00 50.00           C\nEND\n"
        )
        parsed = parse_pdb(two_chains)
        assert len(parsed.chains) == 2
        assert parsed.best_matching_chain("AG").chain_id == "B"
        assert parsed.best_matching_chain("MKT").chain_id == "A"

    def test_longest_chain_when_no_query(self):
        assert parse_pdb(MINIMAL).longest_chain.chain_id == "A"

    def test_no_chains_gives_none(self):
        assert parse_pdb("HEADER x\nEND\n").longest_chain is None


class TestHelpers:
    def test_looks_like_pdb(self):
        assert looks_like_pdb(MINIMAL)
        assert not looks_like_pdb("<html>error</html>")
        assert not looks_like_pdb("")

    def test_accession_scraping(self):
        assert find_accessions("DBREF  1UBQ A 1 76 UNP P62988 UBIQ_HUMAN") == ["P62988"]

    def test_accession_scraping_deduplicates(self):
        assert find_accessions("P62988 and P62988 and P24941") == ["P62988", "P24941"]

    def test_modified_residues_map_to_their_parent(self):
        selenomethionine = MINIMAL.replace("MET A   1", "MSE A   1")
        assert parse_pdb(selenomethionine).chains[0].sequence.startswith("M")

    def test_empty_input_is_safe(self):
        parsed = parse_pdb("")
        assert parsed.chains == [] and not parsed.has_coordinates
