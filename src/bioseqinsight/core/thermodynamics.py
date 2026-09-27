"""Melting-temperature calculation.

Two methods are implemented, and the method actually used is always reported
in the result rather than left implicit:

``wallace``
    :math:`T_m = 2(A+T) + 4(G+C)`. Valid only for short oligonucleotides
    (roughly 14-20 nt). Applying it to gene-length DNA produces numbers in
    the thousands, which v1.0 could still display.

``nearest_neighbor``
    The unified nearest-neighbour parameter set of SantaLucia (1998) with
    helix-initiation terms, terminal A/T penalties, a symmetry correction and
    the SantaLucia salt correction for monovalent cation concentration.

Reference
---------
SantaLucia J. (1998) A unified view of polymer, dumbbell, and oligonucleotide
DNA nearest-neighbor thermodynamics. PNAS 95(4):1460-1465.
"""

from __future__ import annotations

import math

from ..models.results import TmResult
from .alphabet import canonical_dna
from .codon import reverse_complement

# (delta H kcal/mol, delta S cal/(mol*K)) for each nearest-neighbour doublet.
NN_PARAMS: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}

# Helix initiation with a terminal G/C or terminal A/T base pair.
INIT_GC = (0.1, -2.8)
INIT_AT = (2.3, 4.1)
SYMMETRY_CORRECTION_DS = -1.4  # cal/(mol*K), self-complementary duplexes

R = 1.98720425864083  # gas constant, cal/(mol*K)
WALLACE_MAX_LENGTH = 20
NN_MAX_LENGTH = 100  # above this the two-state model stops being appropriate


class ThermodynamicsError(ValueError):
    """Raised when a Tm cannot be computed from the given input."""


def wallace_tm(sequence: str) -> float:
    """Wallace rule value in degrees Celsius."""
    seq = canonical_dna(sequence)
    if not seq:
        raise ThermodynamicsError("No A/C/G/T bases found.")
    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")
    return float(2 * at + 4 * gc)


def is_self_complementary(sequence: str) -> bool:
    seq = canonical_dna(sequence)
    return bool(seq) and seq == reverse_complement(seq)


def nearest_neighbor_thermodynamics(sequence: str) -> tuple[float, float]:
    """Return (delta H in kcal/mol, delta S in cal/(mol*K)) for the duplex."""
    seq = canonical_dna(sequence)
    if len(seq) < 2:
        raise ThermodynamicsError(
            "Nearest-neighbour thermodynamics needs at least 2 canonical bases."
        )
    dh = 0.0
    ds = 0.0
    for i in range(len(seq) - 1):
        doublet = seq[i : i + 2]
        try:
            d_h, d_s = NN_PARAMS[doublet]
        except KeyError as exc:  # pragma: no cover - canonical_dna guarantees ACGT
            raise ThermodynamicsError(f"Unsupported doublet {doublet!r}") from exc
        dh += d_h
        ds += d_s

    for terminal in (seq[0], seq[-1]):
        init_h, init_s = INIT_GC if terminal in "GC" else INIT_AT
        dh += init_h
        ds += init_s

    if is_self_complementary(seq):
        ds += SYMMETRY_CORRECTION_DS
    return dh, ds


def salt_correction(ds: float, length: int, na_mM: float) -> float:
    """SantaLucia (1998) monovalent-cation correction applied to delta S."""
    if na_mM <= 0:
        raise ThermodynamicsError("Sodium concentration must be positive.")
    return ds + 0.368 * (length - 1) * math.log(na_mM / 1000.0)


def nearest_neighbor_tm(
    sequence: str,
    primer_nM: float = 250.0,
    na_mM: float = 50.0,
) -> float:
    """Nearest-neighbour Tm in degrees Celsius.

    ``primer_nM`` is the total strand concentration. For a non
    self-complementary duplex with both strands at equal concentration the
    effective term is ``CT/4``; for a self-complementary duplex it is ``CT``.
    """
    seq = canonical_dna(sequence)
    dh, ds = nearest_neighbor_thermodynamics(seq)
    ds = salt_correction(ds, len(seq), na_mM)
    ct = primer_nM * 1e-9
    if ct <= 0:
        raise ThermodynamicsError("Primer concentration must be positive.")
    effective = ct if is_self_complementary(seq) else ct / 4.0
    denominator = ds + R * math.log(effective)
    if denominator == 0:  # pragma: no cover - numerically unreachable
        raise ThermodynamicsError("Degenerate thermodynamic denominator.")
    return (dh * 1000.0) / denominator - 273.15


def melting_temperature(
    sequence: str,
    primer_nM: float = 250.0,
    na_mM: float = 50.0,
    method: str = "auto",
) -> TmResult:
    """Compute Tm, choosing the appropriate method and saying which was used.

    ``method`` is ``"auto"`` (Wallace at or below 20 nt, nearest-neighbour
    above), ``"wallace"`` or ``"nearest_neighbor"``.
    """
    seq = canonical_dna(sequence)
    if not seq:
        raise ThermodynamicsError(
            "No canonical A/C/G/T bases found; Tm is undefined for this input."
        )
    if method not in ("auto", "wallace", "nearest_neighbor"):
        raise ThermodynamicsError(f"Unknown Tm method {method!r}")

    wallace = wallace_tm(seq)
    chosen = method
    if method == "auto":
        chosen = "wallace" if len(seq) <= WALLACE_MAX_LENGTH else "nearest_neighbor"

    if chosen == "wallace":
        note = None
        if len(seq) > WALLACE_MAX_LENGTH:
            note = (
                f"The Wallace rule was requested for a {len(seq)} nt sequence. "
                "It is only meaningful for roughly 14-20 nt oligonucleotides."
            )
        return TmResult(
            tm_c=wallace,
            method="wallace",
            length=len(seq),
            wallace_c=wallace,
            note=note,
        )

    tm = nearest_neighbor_tm(seq, primer_nM=primer_nM, na_mM=na_mM)
    note = None
    if len(seq) > NN_MAX_LENGTH:
        note = (
            f"Two-state nearest-neighbour thermodynamics was applied to a {len(seq)} nt "
            "sequence. Long duplexes do not melt in a two-state fashion, so treat this "
            "value as an ordering statistic rather than a measurement."
        )
    return TmResult(
        tm_c=tm,
        method="nearest_neighbor",
        length=len(seq),
        wallace_c=wallace,
        na_mM=na_mM,
        primer_nM=primer_nM,
        note=note,
    )


__all__ = [
    "NN_PARAMS",
    "ThermodynamicsError",
    "is_self_complementary",
    "melting_temperature",
    "nearest_neighbor_thermodynamics",
    "nearest_neighbor_tm",
    "salt_correction",
    "wallace_tm",
]
