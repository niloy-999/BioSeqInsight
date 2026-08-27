import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from structure_predictions import (
    _format_ca_score,
    _mean_plddt_from_pdb,
    calculate_hydrophobicity,
    calculate_molecular_weight,
    predict_secondary_structure,
)


UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)


def test_reject_dna_as_protein():
    with pytest.raises(ValueError):
        calculate_hydrophobicity("ATGCATGCATGCATGC")


def test_ubiquitin_mass_and_gravy():
    mass = calculate_molecular_weight(UBIQUITIN)
    gravy = calculate_hydrophobicity(UBIQUITIN)
    assert 8500 < mass < 8700
    assert -0.6 < gravy < 0.2


def test_ss_length_matches():
    data = predict_secondary_structure(UBIQUITIN)
    assert data["length"] == 76
    assert len(data["visual"]) == 76


def test_plddt_scale_0_100():
    pdb = (
        "ATOM      1  CA  MET A   1       0.000   0.000   0.000  1.00 92.50           C\n"
        "ATOM      2  CA  GLN A   2       1.000   0.000   0.000  1.00 80.00           C\n"
    )
    assert abs(_mean_plddt_from_pdb(pdb) - 86.25) < 0.01


def test_score_label_distinguishes_prediction_and_crystal():
    pred = (
        "TITLE     ESMFOLD V1 PREDICTION\n"
        "ATOM      1  CA  MET A   1       0.000   0.000   0.000  1.00 90.00           C\n"
    )
    xtal = (
        "HEADER    KINASE                                  01-JAN-99   1AQ1\n"
        "ATOM      1  CA  MET A   1       0.000   0.000   0.000  1.00 25.00           C\n"
    )
    assert "pLDDT" in _format_ca_score(pred)
    assert "B-factor" in _format_ca_score(xtal)
    assert "not pLDDT" in _format_ca_score(xtal)


def test_plddt_scale_0_1():
    pdb = (
        "ATOM      1  CA  MET A   1       0.000   0.000   0.000  1.00  0.90           C\n"
    )
    assert abs(_mean_plddt_from_pdb(pdb) - 90.0) < 0.01
