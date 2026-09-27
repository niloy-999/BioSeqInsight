"""NCBI genetic code tables, implemented without external dependencies.

Tables 1 (standard), 2 (vertebrate mitochondrial), 4 (mould/protozoan
mitochondrial) and 11 (bacterial/archaeal/plant plastid) are provided; these
cover the organisms represented in the BioSeqInsight benchmark set. The
tables are stored as 64-character strings in the canonical NCBI base order
(TCAG x TCAG x TCAG), which makes them directly checkable against
https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
"""

from __future__ import annotations

from dataclasses import dataclass

_BASES = "TCAG"

# Codons in NCBI order: base1 slowest, base3 fastest.
CODON_ORDER: tuple[str, ...] = tuple(
    b1 + b2 + b3 for b1 in _BASES for b2 in _BASES for b3 in _BASES
)

_NCBI_TABLES: dict[int, tuple[str, str, tuple[str, ...]]] = {
    # id: (name, 64 amino acids in NCBI codon order, start codons)
    1: (
        "Standard",
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        ("TTG", "CTG", "ATG"),
    ),
    2: (
        "Vertebrate Mitochondrial",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        ("ATT", "ATC", "ATA", "ATG", "GTG"),
    ),
    4: (
        "Mold, Protozoan, and Coelenterate Mitochondrial",
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        ("TTA", "TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"),
    ),
    11: (
        "Bacterial, Archaeal and Plant Plastid",
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        ("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"),
    ),
}

STOP = "*"


@dataclass(frozen=True)
class GeneticCode:
    """A single NCBI translation table."""

    table_id: int
    name: str
    forward: dict[str, str]
    starts: frozenset[str]
    stops: frozenset[str]

    def translate_codon(self, codon: str, unknown: str = "X") -> str:
        return self.forward.get(codon.upper(), unknown)

    def is_start(self, codon: str) -> bool:
        return codon.upper() in self.starts

    def is_stop(self, codon: str) -> bool:
        return codon.upper() in self.stops


def _build(table_id: int) -> GeneticCode:
    try:
        name, aas, start_codons = _NCBI_TABLES[table_id]
    except KeyError as exc:  # pragma: no cover - defensive
        raise KeyError(
            f"Genetic code table {table_id} is not bundled. "
            f"Available: {sorted(_NCBI_TABLES)}"
        ) from exc
    if len(aas) != 64:  # pragma: no cover - guards a typo in the table above
        raise ValueError(f"Table {table_id} must define 64 codons, got {len(aas)}")
    forward = dict(zip(CODON_ORDER, aas))
    stop_codons = {codon for codon, aa in forward.items() if aa == STOP}
    unknown_starts = set(start_codons) - set(CODON_ORDER)
    if unknown_starts:  # pragma: no cover - guards a typo in the table above
        raise ValueError(f"Unknown start codons in table {table_id}: {unknown_starts}")
    return GeneticCode(
        table_id=table_id,
        name=name,
        forward=forward,
        starts=frozenset(start_codons),
        stops=frozenset(stop_codons),
    )


_CACHE: dict[int, GeneticCode] = {}


def get_genetic_code(table_id: int = 1) -> GeneticCode:
    """Return a cached :class:`GeneticCode` for an NCBI table id."""
    if table_id not in _CACHE:
        _CACHE[table_id] = _build(table_id)
    return _CACHE[table_id]


def available_tables() -> dict[int, str]:
    return {tid: value[0] for tid, value in _NCBI_TABLES.items()}


COMPLEMENT = {
    "A": "T", "T": "A", "G": "C", "C": "G", "U": "A", "N": "N",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D",
}

_RC_TABLE = str.maketrans(
    "".join(COMPLEMENT) + "".join(COMPLEMENT).lower(),
    "".join(COMPLEMENT.values()) + "".join(COMPLEMENT.values()).lower(),
)


def complement(sequence: str) -> str:
    """Complement a nucleotide string, preserving IUPAC ambiguity codes."""
    return sequence.translate(_RC_TABLE)


def reverse_complement(sequence: str) -> str:
    """Reverse complement, preserving IUPAC ambiguity codes."""
    return complement(sequence)[::-1]


__all__ = [
    "CODON_ORDER",
    "STOP",
    "GeneticCode",
    "available_tables",
    "complement",
    "get_genetic_code",
    "reverse_complement",
]
