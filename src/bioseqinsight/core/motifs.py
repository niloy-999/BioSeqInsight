"""Motif searching.

This is exact pattern matching, not statistical motif discovery: BioSeqInsight
is not MEME and does not claim to be. Three query styles are supported and the
style used is recorded in every hit so that exported tables are unambiguous:

* a literal string (``GAATTC``);
* an IUPAC degenerate string (``GGNTCA``, ``WGATAR``);
* a raw regular expression, when ``regex=True``.

Overlapping occurrences are reported, which plain ``str.find`` loops in v1.0
did not guarantee for degenerate patterns.
"""

from __future__ import annotations

import re

from ..models.results import MotifHit
from .alphabet import DNA_BASES, clean_dna
from .codon import reverse_complement

IUPAC_EXPANSION: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}


class MotifError(ValueError):
    """Raised when a motif pattern cannot be compiled."""


def iupac_to_regex(motif: str) -> str:
    """Translate an IUPAC degenerate motif into a regular expression."""
    pattern = clean_dna(motif)
    if not pattern:
        raise MotifError("Empty motif: enter at least one nucleotide code.")
    parts: list[str] = []
    for ch in pattern:
        expansion = IUPAC_EXPANSION.get(ch)
        if expansion is None:
            raise MotifError(f"{ch!r} is not an IUPAC nucleotide code.")
        parts.append(expansion if len(expansion) == 1 else f"[{expansion}]")
    return "".join(parts)


def compile_motif(motif: str, regex: bool = False) -> re.Pattern[str]:
    source = motif if regex else iupac_to_regex(motif)
    try:
        return re.compile(source, re.IGNORECASE)
    except re.error as exc:
        raise MotifError(f"Invalid pattern {motif!r}: {exc}") from exc


def find_motif(
    sequence: str,
    motif: str,
    both_strands: bool = False,
    regex: bool = False,
) -> list[MotifHit]:
    """Return every (overlapping) occurrence of ``motif``.

    With ``both_strands=True`` the reverse strand is also scanned and hits are
    reported with coordinates on the *forward* strand so the two sets can be
    compared directly.
    """
    seq = clean_dna(sequence)
    if not seq:
        return []
    pattern = compile_motif(motif, regex=regex)
    label = motif.upper() if not regex else motif

    hits = [
        MotifHit(motif=label, start=m.start(), end=m.start() + _span(m), strand="+")
        for m in _overlapping(pattern, seq)
    ]

    if both_strands:
        rev = reverse_complement(seq)
        n = len(seq)
        for m in _overlapping(pattern, rev):
            width = _span(m)
            start = n - (m.start() + width)
            hits.append(MotifHit(motif=label, start=start, end=start + width, strand="-"))
        hits.sort(key=lambda h: (h.start, h.strand))
    return hits


def _overlapping(pattern: re.Pattern[str], seq: str):
    pos = 0
    while pos <= len(seq):
        match = pattern.search(seq, pos)
        if match is None:
            return
        yield match
        pos = match.start() + 1


def _span(match: re.Match[str]) -> int:
    return max(match.end() - match.start(), 1)


def gc_skew(sequence: str, window: int = 100) -> list[tuple[int, float]]:
    """Sliding-window GC skew, ``(G - C) / (G + C)``.

    Included because replication-origin teaching examples are a common use of
    the sequence tab; the window is returned with its start coordinate so the
    series can be plotted or exported directly.
    """
    seq = clean_dna(sequence)
    if window < 1:
        raise ValueError("window must be >= 1")
    out: list[tuple[int, float]] = []
    for start in range(0, max(len(seq) - window + 1, 0), window):
        chunk = seq[start : start + window]
        g, c = chunk.count("G"), chunk.count("C")
        out.append((start, (g - c) / (g + c) if (g + c) else 0.0))
    return out


def restriction_summary(sequence: str, sites: dict[str, str]) -> dict[str, list[int]]:
    """Count user-supplied recognition sites.

    ``sites`` maps an enzyme name to an IUPAC recognition sequence. No enzyme
    database is bundled: BioSeqInsight does not redistribute REBASE.
    """
    seq = clean_dna(sequence)
    return {
        name: [h.start for h in find_motif(seq, site)]
        for name, site in sites.items()
        if set(site.upper()) <= set(IUPAC_EXPANSION)
    }


__all__ = [
    "DNA_BASES",
    "IUPAC_EXPANSION",
    "MotifError",
    "compile_motif",
    "find_motif",
    "gc_skew",
    "iupac_to_regex",
    "restriction_summary",
]
