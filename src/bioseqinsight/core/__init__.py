"""Deterministic, dependency-free computational core."""

from .alphabet import (
    FastaRecord,
    SequenceError,
    clean_dna,
    clean_protein,
    detect_alphabet,
    parse_fasta,
    write_fasta,
)
from .protein import ProteinAnalyzer
from .sequence import SequenceAnalyzer, analyze_sequence, gc_content

__all__ = [
    "FastaRecord",
    "ProteinAnalyzer",
    "SequenceAnalyzer",
    "SequenceError",
    "analyze_sequence",
    "clean_dna",
    "clean_protein",
    "detect_alphabet",
    "gc_content",
    "parse_fasta",
    "write_fasta",
]
