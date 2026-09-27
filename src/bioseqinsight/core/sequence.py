"""The nucleotide analysis facade used by every front end.

``SequenceAnalyzer`` is the clean software API that the v2.0 architecture puts
between the GUI and the computation. The GUI calls::

    result = SequenceAnalyzer(text, identifier="seq1").analyze()

and formats ``result``; it never computes GC content, translates DNA or
handles errors itself.
"""

from __future__ import annotations

from .. import __version__
from ..models.results import ORF, MotifHit, SequenceResult, TmResult
from .alphabet import (
    DNA_BASES,
    canonical_dna,
    clean_dna,
    validate_dna,
)
from .codon import reverse_complement as _rc
from .motifs import find_motif, gc_skew
from .orf import find_orfs
from .thermodynamics import melting_temperature
from .translation import six_frame_translation, transcribe, translate_peptides


def gc_content(sequence: str) -> float:
    """GC percentage over canonical bases only.

    Ambiguity codes are excluded from both numerator and denominator so that
    a sequence padded with N does not silently depress the reported GC%.
    """
    seq = canonical_dna(sequence)
    if not seq:
        return 0.0
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def nucleotide_composition(sequence: str) -> dict[str, int]:
    """Counts of A/C/G/T plus a single ``N`` bucket for ambiguity codes."""
    seq = clean_dna(sequence)
    counts = {base: seq.count(base) for base in "ATGC"}
    counts["N"] = len(seq) - sum(counts.values())
    return counts


class SequenceAnalyzer:
    """Deterministic analysis of one nucleotide sequence."""

    def __init__(self, sequence: str, identifier: str = "sequence", table: int = 1):
        self.identifier = identifier
        self.table = table
        self.raw = sequence
        self.sequence = validate_dna(sequence)
        self.canonical = canonical_dna(self.sequence)

    def __len__(self) -> int:
        return len(self.sequence)

    # -- individual operations ---------------------------------------------

    def get_length(self) -> int:
        return len(self.sequence)

    def get_composition(self) -> dict[str, int]:
        return nucleotide_composition(self.sequence)

    def get_gc_content(self) -> float:
        return gc_content(self.sequence)

    def get_gc_skew(self, window: int = 100) -> list[tuple[int, float]]:
        return gc_skew(self.sequence, window=window)

    def get_reverse_complement(self) -> str:
        return _rc(self.sequence)

    def transcribe(self) -> str:
        return transcribe(self.sequence)

    def translate(self, frame: int = 0) -> list[str]:
        return translate_peptides(self.sequence, table=self.table, frame=frame)

    def six_frames(self) -> dict[str, str]:
        return six_frame_translation(self.sequence, table=self.table)

    def find_motifs(self, motif: str, both_strands: bool = True, regex: bool = False) -> list[MotifHit]:
        return find_motif(self.sequence, motif, both_strands=both_strands, regex=regex)

    def find_orfs(self, min_aa: int = 30, include_reverse: bool = True) -> list[ORF]:
        return find_orfs(
            self.sequence,
            min_aa=min_aa,
            table=self.table,
            include_reverse=include_reverse,
        )

    def calculate_tm(self, primer_nM: float = 250.0, na_mM: float = 50.0) -> TmResult:
        return melting_temperature(self.sequence, primer_nM=primer_nM, na_mM=na_mM)

    # -- full report --------------------------------------------------------

    def analyze(
        self,
        motif: str | None = None,
        min_orf_aa: int = 30,
        include_reverse_orfs: bool = True,
        include_translation: bool = True,
        include_sequences: bool = True,
    ) -> SequenceResult:
        """Run the standard analysis battery and return a structured result.

        ``include_sequences=False`` omits the reverse complement, mRNA and
        peptide strings. Batch runs over thousands of sequences use that to
        keep exported tables small.
        """
        composition = self.get_composition()
        n = len(self.sequence)
        warnings: list[str] = []
        ambiguous = composition["N"]
        if ambiguous:
            warnings.append(
                f"{ambiguous} ambiguity code(s) present; GC%, Tm and translation "
                "use canonical A/C/G/T positions only."
            )
        if n < 3:
            warnings.append("Sequence is shorter than one codon; translation is empty.")

        tm: TmResult | None = None
        if self.canonical:
            tm = self.calculate_tm()

        return SequenceResult(
            identifier=self.identifier,
            length=n,
            alphabet="dna",
            gc_percent=self.get_gc_content(),
            composition=composition,
            composition_percent={
                base: (100.0 * count / n if n else 0.0)
                for base, count in composition.items()
            },
            ambiguous_count=ambiguous,
            reverse_complement=self.get_reverse_complement() if include_sequences else None,
            mrna=self.transcribe() if include_sequences else None,
            peptides=(self.translate() if include_translation and include_sequences else []),
            orfs=self.find_orfs(min_aa=min_orf_aa, include_reverse=include_reverse_orfs),
            motifs=(self.find_motifs(motif) if motif else []),
            tm=tm,
            warnings=warnings,
            software_version=__version__,
        )


def analyze_sequence(sequence: str, identifier: str = "sequence", **kwargs) -> SequenceResult:
    """Functional shortcut: ``analyze_sequence(text)`` -> :class:`SequenceResult`."""
    return SequenceAnalyzer(sequence, identifier=identifier).analyze(**kwargs)


__all__ = [
    "DNA_BASES",
    "SequenceAnalyzer",
    "analyze_sequence",
    "gc_content",
    "nucleotide_composition",
]
