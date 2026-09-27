"""Cross-validation of the dependency-free core against Biopython.

BioSeqInsight v2.0 implements its sequence and protein calculations in plain
Python so that the whole package can be installed and tested anywhere,
including on a machine with no compiler and no network. That choice is only
defensible if the implementations agree with an established reference, so
this module checks them against Biopython wherever Biopython is installed.

These tests are skipped, not failed, when Biopython is absent: it is an
optional dependency. In continuous integration at least one matrix job
installs it, so the comparison always runs somewhere before a release.

Where an exact match is not expected the tolerance is stated and the reason
is documented in the test itself.
"""

from __future__ import annotations

import pytest

Bio = pytest.importorskip("Bio")

from Bio.Seq import Seq  # noqa: E402
from Bio.SeqUtils import gc_fraction  # noqa: E402
from Bio.SeqUtils.ProtParam import ProteinAnalysis  # noqa: E402

from bioseqinsight.core.codon import reverse_complement  # noqa: E402
from bioseqinsight.core.protein import (  # noqa: E402
    aromaticity,
    gravy,
    isoelectric_point,
    molecular_weight,
)
from bioseqinsight.core.sequence import gc_content  # noqa: E402
from bioseqinsight.core.translation import translate  # noqa: E402

DNA_CASES = [
    "ATGGCTAGCTAGGCTTACCGAT",
    "GGGGCCCCAAAATTTT",
    "ATGAAACCCGGGTTTTAA",
    "ACGTACGTACGTACGTACGTACGTACGT",
]

PROTEIN_CASES = [
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ",
    "ACDEFGHIKLMNPQRSTVWY",
    "GIVEQCCTSICSLYQLENYCN",
]


class TestNucleotideParity:
    @pytest.mark.parametrize("sequence", DNA_CASES)
    def test_reverse_complement_matches(self, sequence):
        assert reverse_complement(sequence) == str(Seq(sequence).reverse_complement())

    @pytest.mark.parametrize("sequence", DNA_CASES)
    def test_gc_content_matches(self, sequence):
        assert gc_content(sequence) == pytest.approx(100.0 * gc_fraction(sequence), abs=1e-6)

    @pytest.mark.parametrize("sequence", DNA_CASES)
    def test_translation_matches(self, sequence):
        usable = sequence[: len(sequence) - len(sequence) % 3]
        assert translate(usable) == str(Seq(usable).translate())

    def test_bacterial_table_translation_matches(self):
        sequence = "ATGAAACCCGGGTGATAA"
        assert translate(sequence, table=11) == str(Seq(sequence).translate(table=11))


class TestProteinParity:
    @pytest.mark.parametrize("sequence", PROTEIN_CASES)
    def test_molecular_weight_matches(self, sequence):
        # Both use the same IUPAC average residue masses, so agreement should
        # be to within floating-point noise.
        assert molecular_weight(sequence) == pytest.approx(
            ProteinAnalysis(sequence).molecular_weight(), abs=1e-4
        )

    @pytest.mark.parametrize("sequence", PROTEIN_CASES)
    def test_gravy_matches(self, sequence):
        assert gravy(sequence) == pytest.approx(ProteinAnalysis(sequence).gravy(), abs=1e-6)

    @pytest.mark.parametrize("sequence", PROTEIN_CASES)
    def test_aromaticity_matches(self, sequence):
        assert aromaticity(sequence) == pytest.approx(
            ProteinAnalysis(sequence).aromaticity(), abs=1e-6
        )




class TestDocumentedDifferences:
    def test_melting_temperature_is_close_but_not_identical(self):
        """Tm agrees with Biopython to within a few degrees, by design.

        Both implement SantaLucia (1998) nearest-neighbour thermodynamics,
        but Biopython's ``Tm_NN`` defaults to the Owczarzy (2004) salt
        correction and 25 nM strand concentration, whereas BioSeqInsight uses
        the SantaLucia (1998) correction and 250 nM. The conditions actually
        used are reported in every ``TmResult``, so the two are comparable
        once matched. This test documents the expected size of the gap rather
        than asserting equality.
        """
        from Bio.SeqUtils import MeltingTemp as mt

        from bioseqinsight.core.thermodynamics import nearest_neighbor_tm

        for sequence in ("ATGCATGCATGCATGCATGC", "GCGCGCGCATATATATGCGC"):
            ours = nearest_neighbor_tm(sequence, primer_nM=250.0, na_mM=50.0)
            theirs = float(mt.Tm_NN(Seq(sequence)))
            assert abs(ours - theirs) < 10.0

    def test_ambiguity_handling_differs_deliberately(self):
        """We translate ambiguous codons to X; Biopython may raise or guess.

        Keeping the residue in place preserves the reading frame, which is the
        behaviour a user expects when a sequence contains an N.
        """
        assert translate("ATGNNNAAA") == "MXK"

    @pytest.mark.parametrize("sequence", PROTEIN_CASES)
    def test_isoelectric_point_is_close_but_not_identical(self, sequence):
        """pI disagrees with Biopython by up to a few tenths of a pH unit.

        This was not the intended outcome. The docstring in
        ``core/protein.py`` originally claimed BioSeqInsight used "the pKa
        set used by EMBOSS/Biopython", and this test originally asserted
        near-exact agreement (abs=0.05) on that basis. Running it against a
        real Biopython installation showed differences of 0.23-0.26 pH
        units for two of the four test sequences — far larger than the
        bisection tolerance could explain, which means the pKa tables
        themselves differ, not just the numerical method.

        There is no single agreed pKa scale for peptide ionisable groups:
        independent pI calculators are documented to disagree by several
        tenths of a pH unit, sometimes more, depending on which published
        pKa set they use (see Kozlowski, 2016, "IPC - Isoelectric Point
        Calculator", Biology Direct 11:55, which benchmarks this exact
        disagreement across common tools). BioSeqInsight's pKa table is
        printed alongside the result's method where relevant; a reader who
        needs pI to agree with a specific external tool should compare pKa
        tables rather than assume any two calculators agree.

        The tolerance below is wide enough to accept a difference of this
        known size while still catching a genuinely broken calculation (for
        example, a sign error, which would put the two several pH units
        apart rather than a few tenths).
        """
        ours = isoelectric_point(sequence)
        theirs = ProteinAnalysis(sequence).isoelectric_point()
        assert ours == pytest.approx(theirs, abs=0.5)

    def test_isoelectric_point_is_self_consistent(self):
        """Independent of Biopython: charge is (near) zero at our own pI.

        This is the check that actually matters for correctness. Agreement
        with another tool's specific pKa table is a matter of convention;
        the net charge vanishing at the pH the bisection converged on is a
        property the calculation must have regardless of which table is
        used.
        """
        from bioseqinsight.core.protein import charge_at_ph

        for sequence in PROTEIN_CASES:
            pi = isoelectric_point(sequence)
            assert abs(charge_at_ph(sequence, pi)) < 1e-2
