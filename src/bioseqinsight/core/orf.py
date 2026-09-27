"""Open reading frame detection.

The v1.0 implementation reported the first start codon it met and could emit
overlapping duplicates of the same ORF. v2.0 defines the behaviour precisely
and tests it:

* Six frames are scanned by default (three forward, three reverse).
* An ORF runs from a start codon to the first in-frame stop codon.
* Only the *longest* ORF per (strand, frame, stop position) is reported, so
  internal downstream starts do not create nested duplicates.
* ORFs that reach the end of the sequence without a stop codon are reported
  and flagged ``unterminated``.
* Ambiguity codes stay in place; a codon containing one is neither a start
  nor a stop and translates to ``X``. Removing them would shift the frame.
* Coordinates are 0-based half-open on the strand that was scanned; the
  helper :func:`to_forward_coordinates` maps reverse hits back onto the
  forward strand for export.
"""

from __future__ import annotations

from ..models.results import ORF
from .alphabet import DNA_BASES, clean_dna
from .codon import get_genetic_code, reverse_complement


def find_orfs(
    sequence: str,
    min_aa: int = 30,
    table: int = 1,
    include_reverse: bool = True,
    require_atg: bool = True,
    include_unterminated: bool = True,
) -> list[ORF]:
    """Return ORFs sorted by decreasing peptide length.

    ``require_atg`` restricts starts to ATG (the default, and what most users
    expect); setting it to ``False`` uses every alternative start codon
    defined by the selected NCBI table, which matters for bacterial genomes.
    """
    if min_aa < 1:
        raise ValueError("min_aa must be >= 1")
    code = get_genetic_code(table)
    forward = clean_dna(sequence)
    strands: list[tuple[str, str]] = [("+", forward)]
    if include_reverse:
        strands.append(("-", reverse_complement(forward)))

    starts = {"ATG"} if require_atg else set(code.starts)
    results: list[ORF] = []

    for strand, seq in strands:
        for frame in range(3):
            open_start: int | None = None
            for pos in range(frame, len(seq) - 2, 3):
                codon = seq[pos : pos + 3]
                if open_start is None:
                    if codon in starts:
                        open_start = pos
                    continue
                if code.is_stop(codon):
                    aa_len = (pos - open_start) // 3
                    if aa_len >= min_aa:
                        results.append(
                            _build(code, strand, frame, seq, open_start, pos + 3, aa_len, False)
                        )
                    open_start = None
            if open_start is not None and include_unterminated:
                end = len(seq) - ((len(seq) - open_start) % 3)
                aa_len = (end - open_start) // 3
                if aa_len >= min_aa:
                    results.append(
                        _build(code, strand, frame, seq, open_start, end, aa_len, True)
                    )

    results.sort(key=lambda o: (-o.aa_length, o.strand, o.start))
    return results


def _build(code, strand, frame, seq, start, end, aa_len, unterminated) -> ORF:
    coding = seq[start : start + aa_len * 3]
    peptide = "M" + "".join(
        code.translate_codon(coding[i : i + 3]) if set(coding[i : i + 3]) <= DNA_BASES else "X"
        for i in range(3, len(coding), 3)
    )
    return ORF(
        strand=strand,
        frame=frame + 1,
        start=start,
        end=end,
        aa_length=aa_len,
        protein=peptide,
        unterminated=unterminated,
    )


def to_forward_coordinates(orf: ORF, sequence_length: int) -> tuple[int, int]:
    """Map an ORF's coordinates back onto the forward strand.

    Reverse-strand coordinates are measured on the reverse complement, which
    is convenient while scanning but unhelpful in an exported table.
    """
    if orf.strand == "+":
        return orf.start, orf.end
    return sequence_length - orf.end, sequence_length - orf.start


def longest_orf(sequence: str, **kwargs) -> ORF | None:
    """Convenience wrapper returning only the longest ORF, or ``None``."""
    orfs = find_orfs(sequence, **kwargs)
    return orfs[0] if orfs else None


__all__ = ["find_orfs", "longest_orf", "to_forward_coordinates"]
