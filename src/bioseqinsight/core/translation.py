"""Deterministic translation of nucleotide sequences."""

from __future__ import annotations

from .alphabet import DNA_BASES, clean_dna
from .codon import STOP, get_genetic_code, reverse_complement


def translate(
    sequence: str,
    table: int = 1,
    frame: int = 0,
    to_stop: bool = False,
    unknown: str = "X",
) -> str:
    """Translate one reading frame and return a single peptide string.

    ``frame`` is 0-based (0, 1 or 2). Incomplete trailing codons are ignored.
    Stop codons are rendered as ``*`` unless ``to_stop`` truncates at the
    first one.

    Ambiguity codes are **kept in place** and their codons translate to
    ``unknown`` (``X`` by default). Deleting them instead, as earlier versions
    did, silently shifts the reading frame for every downstream codon and
    produces a plausible-looking but wrong peptide.
    """
    if frame not in (0, 1, 2):
        raise ValueError("frame must be 0, 1 or 2")
    code = get_genetic_code(table)
    seq = clean_dna(sequence)[frame:]
    usable = seq[: len(seq) - (len(seq) % 3)]
    out: list[str] = []
    for i in range(0, len(usable), 3):
        codon = usable[i : i + 3]
        aa = (
            code.translate_codon(codon, unknown=unknown)
            if set(codon) <= DNA_BASES
            else unknown
        )
        if aa == STOP and to_stop:
            break
        out.append(aa)
    return "".join(out)


def translate_peptides(sequence: str, table: int = 1, frame: int = 0) -> list[str]:
    """Translate a frame and split the product on stop codons.

    Empty fragments (consecutive stops) are dropped, so the result is the list
    of uninterrupted peptides in that frame.
    """
    protein = translate(sequence, table=table, frame=frame)
    return [part for part in protein.split(STOP) if part]


def six_frame_translation(sequence: str, table: int = 1) -> dict[str, str]:
    """Return all six frames keyed ``"+1".."+3"`` and ``"-1".."-3"``."""
    seq = clean_dna(sequence)
    rev = reverse_complement(seq)
    frames: dict[str, str] = {}
    for f in range(3):
        frames[f"+{f + 1}"] = translate(seq, table=table, frame=f)
        frames[f"-{f + 1}"] = translate(rev, table=table, frame=f)
    return frames


def transcribe(sequence: str) -> str:
    """DNA (coding strand) to mRNA."""
    return clean_dna(sequence).replace("T", "U")


def back_transcribe(sequence: str) -> str:
    """mRNA back to the DNA coding strand."""
    return clean_dna(sequence)


__all__ = [
    "back_transcribe",
    "six_frame_translation",
    "transcribe",
    "translate",
    "translate_peptides",
]
