"""Protein physicochemical analysis.

All descriptors are implemented in plain Python so that they can be tested
offline and so that BioSeqInsight has no hard dependency on a heavyweight
stack. Where Biopython is installed, ``tests/test_biopython_parity.py``
cross-checks molecular weight, isoelectric point, GRAVY and aromaticity
against ``Bio.SeqUtils.ProtParam`` to confirm the implementations agree.

Constants and references
------------------------
* Average residue masses: IUPAC average isotopic masses (the same values
  Biopython uses in ``IUPACData.protein_weights``).
* Hydropathy: Kyte & Doolittle (1982) J Mol Biol 157:105-132.
* Isoelectric point: a Bjellqvist-style iterative charge model. The specific
  pKa values below are one of several sets in common use; independent pI
  calculators are known to disagree by several tenths of a pH unit because
  there is no single agreed scale (Kozlowski, 2016, Biology Direct 11:55).
  ``tests/test_biopython_parity.py`` checks agreement with Biopython within
  that known tolerance, and separately checks that the net charge is zero at
  the reported pI regardless of which table is used — that self-consistency
  check, not agreement with any one external tool, is what correctness means
  here.
* Molar extinction coefficient: Gill & von Hippel (1989) Anal Biochem
  182:319-326.
* Aliphatic index: Ikai (1980) J Biochem 88:1895-1898.
* Secondary-structure propensities: Chou & Fasman (1978) style single-residue
  preferences, used here only as a smoothed sketch (see
  :func:`secondary_structure_sketch`).
"""

from __future__ import annotations

from .. import __version__
from ..models.results import ProteinResult, SecondaryStructureSketch
from .alphabet import PROTEIN_LETTERS, SequenceError, clean_protein, validate_protein

WATER = 18.0153

AVERAGE_RESIDUE_MASS: dict[str, float] = {
    "A": 89.0932, "C": 121.1582, "D": 133.1027, "E": 147.1293, "F": 165.1891,
    "G": 75.0666, "H": 155.1546, "I": 131.1729, "K": 146.1876, "L": 131.1729,
    "M": 149.2113, "N": 132.1179, "P": 115.1305, "Q": 146.1445, "R": 174.2010,
    "S": 105.0926, "T": 119.1192, "V": 117.1463, "W": 204.2252, "Y": 181.1885,
}

KYTE_DOOLITTLE: dict[str, float] = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4, "H": -3.2,
    "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5,
    "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}

POSITIVE_PKS = {"Nterm": 7.5, "K": 10.0, "R": 12.0, "H": 5.98}
NEGATIVE_PKS = {"Cterm": 3.55, "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}

EXTINCTION_TRP = 5500
EXTINCTION_TYR = 1490
EXTINCTION_CYSTINE = 125

HELIX_PROPENSITY: dict[str, float] = {
    "A": 1.42, "C": 0.70, "D": 1.01, "E": 1.51, "F": 1.13, "G": 0.57, "H": 1.00,
    "I": 1.08, "K": 1.16, "L": 1.21, "M": 1.45, "N": 0.67, "P": 0.57, "Q": 1.11,
    "R": 0.98, "S": 0.77, "T": 0.83, "V": 1.06, "W": 1.08, "Y": 0.69,
}
SHEET_PROPENSITY: dict[str, float] = {
    "A": 0.83, "C": 1.19, "D": 0.54, "E": 0.37, "F": 1.38, "G": 0.75, "H": 0.87,
    "I": 1.60, "K": 0.74, "L": 1.30, "M": 1.05, "N": 0.89, "P": 0.55, "Q": 1.10,
    "R": 0.93, "S": 0.75, "T": 1.19, "V": 1.70, "W": 1.37, "Y": 1.47,
}
COIL_PROPENSITY: dict[str, float] = {
    "A": 0.66, "C": 1.19, "D": 1.46, "E": 0.74, "F": 0.60, "G": 1.52, "H": 0.87,
    "I": 0.47, "K": 1.20, "L": 0.59, "M": 0.60, "N": 1.56, "P": 1.52, "Q": 0.92,
    "R": 0.93, "S": 1.43, "T": 1.20, "V": 0.61, "W": 0.60, "Y": 1.14,
}

_SS_GLYPHS = {"H": "H", "E": "E", "C": "-"}


def _standard_only(sequence: str) -> str:
    return "".join(ch for ch in sequence if ch in PROTEIN_LETTERS)


def molecular_weight(sequence: str) -> float:
    """Average isotopic molecular weight in daltons.

    Ambiguity codes (B, J, O, U, X, Z) carry no defined average mass and are
    excluded; :func:`analyze` reports how many residues were skipped so the
    omission is never silent.
    """
    seq = _standard_only(validate_protein(sequence))
    if not seq:
        raise SequenceError("Sequence contains no standard amino acids.")
    return sum(AVERAGE_RESIDUE_MASS[aa] for aa in seq) - (len(seq) - 1) * WATER


def gravy(sequence: str) -> float:
    """Grand average of hydropathy (mean Kyte-Doolittle value)."""
    seq = _standard_only(validate_protein(sequence))
    if not seq:
        raise SequenceError("Sequence contains no standard amino acids.")
    return sum(KYTE_DOOLITTLE[aa] for aa in seq) / len(seq)


def hydropathy_profile(sequence: str, window: int = 9) -> list[float]:
    """Sliding-window Kyte-Doolittle profile (one value per window position)."""
    if window < 1:
        raise ValueError("window must be >= 1")
    seq = _standard_only(validate_protein(sequence))
    if len(seq) < window:
        return []
    values = [KYTE_DOOLITTLE[aa] for aa in seq]
    return [
        sum(values[i : i + window]) / window for i in range(len(values) - window + 1)
    ]


def amino_acid_composition(sequence: str) -> dict[str, int]:
    seq = clean_protein(sequence)
    return {aa: seq.count(aa) for aa in sorted(PROTEIN_LETTERS)}


def aromaticity(sequence: str) -> float:
    """Relative frequency of F, W and Y (Lobry & Gautier, 1994)."""
    seq = _standard_only(validate_protein(sequence))
    if not seq:
        return 0.0
    return sum(seq.count(aa) for aa in "FWY") / len(seq)


def aliphatic_index(sequence: str) -> float:
    """Ikai (1980) aliphatic index, a proxy for thermostability."""
    seq = _standard_only(validate_protein(sequence))
    if not seq:
        return 0.0
    n = len(seq)
    a = 100.0 * seq.count("A") / n
    v = 100.0 * seq.count("V") / n
    il = 100.0 * (seq.count("I") + seq.count("L")) / n
    return a + 2.9 * v + 3.9 * il


def extinction_coefficient(sequence: str) -> tuple[int, int]:
    """Molar extinction at 280 nm: (all Cys reduced, all Cys as cystines)."""
    seq = _standard_only(validate_protein(sequence))
    base = EXTINCTION_TRP * seq.count("W") + EXTINCTION_TYR * seq.count("Y")
    cystines = seq.count("C") // 2
    return base, base + EXTINCTION_CYSTINE * cystines


def charge_at_ph(sequence: str, ph: float) -> float:
    """Net charge of the peptide at a given pH."""
    seq = _standard_only(validate_protein(sequence))
    if not seq:
        return 0.0
    positive = 0.0
    for group, pk in POSITIVE_PKS.items():
        count = 1 if group == "Nterm" else seq.count(group)
        if count:
            positive += count * (1.0 / (1.0 + 10.0 ** (ph - pk)))
    negative = 0.0
    for group, pk in NEGATIVE_PKS.items():
        count = 1 if group == "Cterm" else seq.count(group)
        if count:
            negative += count * (1.0 / (1.0 + 10.0 ** (pk - ph)))
    return positive - negative


def isoelectric_point(sequence: str, tolerance: float = 1e-4) -> float:
    """pI found by bisection on the net-charge curve between pH 0 and 14."""
    seq = _standard_only(validate_protein(sequence))
    if not seq:
        raise SequenceError("Sequence contains no standard amino acids.")
    low, high = 0.0, 14.0
    for _ in range(200):
        mid = (low + high) / 2.0
        charge = charge_at_ph(seq, mid)
        if abs(charge) < tolerance or (high - low) < tolerance:
            return mid
        if charge > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0  # pragma: no cover - bisection always converges


def secondary_structure_sketch(sequence: str, window: int = 5) -> SecondaryStructureSketch:
    """Smoothed single-residue propensity sketch.

    This is explicitly **not** a secondary-structure prediction. It assigns
    each residue the state with the highest Chou-Fasman-style single-residue
    preference and then applies a majority filter over ``window`` residues.
    It is provided for teaching illustration only and the software labels it
    as such everywhere it is displayed or exported.
    """
    seq = _standard_only(validate_protein(sequence))
    if not seq:
        raise SequenceError("Sequence contains no standard amino acids.")
    raw = [
        max(
            ("H", "E", "C"),
            key=lambda state: {
                "H": HELIX_PROPENSITY[aa],
                "E": SHEET_PROPENSITY[aa],
                "C": COIL_PROPENSITY[aa],
            }[state],
        )
        for aa in seq
    ]
    half = max(window // 2, 0)
    smoothed: list[str] = []
    for i in range(len(raw)):
        chunk = raw[max(0, i - half) : i + half + 1]
        smoothed.append(max(("H", "E", "C"), key=chunk.count))
    n = len(smoothed)
    return SecondaryStructureSketch(
        length=n,
        visual="".join(_SS_GLYPHS[s] for s in smoothed),
        helix_percent=100.0 * smoothed.count("H") / n,
        sheet_percent=100.0 * smoothed.count("E") / n,
        coil_percent=100.0 * smoothed.count("C") / n,
        method=(
            f"smoothed single-residue propensities, window={window} "
            "(illustrative sketch, not a structure prediction)"
        ),
    )


class ProteinAnalyzer:
    """Object-oriented facade over the protein functions.

    Every method is independently callable and returns structured data, which
    is what makes the GUI, the CLI and the benchmark share one code path.
    """

    def __init__(self, sequence: str, identifier: str = "protein"):
        self.identifier = identifier
        self.raw = sequence
        self.sequence = validate_protein(sequence)
        self.standard = _standard_only(self.sequence)

    def __len__(self) -> int:
        return len(self.sequence)

    def get_length(self) -> int:
        return len(self.sequence)

    def get_molecular_weight(self) -> float:
        return molecular_weight(self.sequence)

    def get_amino_acid_composition(self) -> dict[str, int]:
        return amino_acid_composition(self.sequence)

    def get_gravy(self) -> float:
        return gravy(self.sequence)

    def get_hydropathy_profile(self, window: int = 9) -> list[float]:
        return hydropathy_profile(self.sequence, window=window)

    def get_isoelectric_point(self) -> float:
        return isoelectric_point(self.sequence)

    def get_secondary_structure_propensity(self, window: int = 5) -> SecondaryStructureSketch:
        return secondary_structure_sketch(self.sequence, window=window)

    def analyze(self, include_sketch: bool = True) -> ProteinResult:
        """Run the full descriptor set and return a :class:`ProteinResult`."""
        composition = amino_acid_composition(self.sequence)
        n = len(self.sequence)
        unknown = n - len(self.standard)
        warnings: list[str] = []
        if unknown:
            warnings.append(
                f"{unknown} residue(s) use ambiguity codes and were excluded from "
                "mass, charge and hydropathy calculations."
            )
        reduced, cystines = extinction_coefficient(self.sequence)
        return ProteinResult(
            identifier=self.identifier,
            length=n,
            molecular_weight=molecular_weight(self.sequence),
            gravy=gravy(self.sequence),
            aromaticity=aromaticity(self.sequence),
            isoelectric_point=isoelectric_point(self.sequence),
            aliphatic_index=aliphatic_index(self.sequence),
            extinction_coefficient_reduced=reduced,
            extinction_coefficient_cystines=cystines,
            composition=composition,
            composition_percent={
                aa: (100.0 * count / n if n else 0.0) for aa, count in composition.items()
            },
            unknown_residues=unknown,
            secondary_structure=(
                secondary_structure_sketch(self.sequence) if include_sketch else None
            ),
            warnings=warnings,
            software_version=__version__,
        )


__all__ = [
    "AVERAGE_RESIDUE_MASS",
    "KYTE_DOOLITTLE",
    "ProteinAnalyzer",
    "aliphatic_index",
    "amino_acid_composition",
    "aromaticity",
    "charge_at_ph",
    "extinction_coefficient",
    "gravy",
    "hydropathy_profile",
    "isoelectric_point",
    "molecular_weight",
    "secondary_structure_sketch",
]
