"""Input normalisation, alphabet validation and FASTA handling.

Every public entry point in BioSeqInsight funnels user input through this
module first. Normalisation is therefore defined in exactly one place and is
covered by its own unit tests, which removes a whole class of
"works-in-the-GUI-only" behaviour present in v1.0.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

DNA_BASES = frozenset("ACGT")
RNA_BASES = frozenset("ACGU")
# IUPAC nucleotide ambiguity codes.
DNA_AMBIGUITY = frozenset("RYSWKMBDHVN")
DNA_ALPHABET = DNA_BASES | DNA_AMBIGUITY
PROTEIN_LETTERS = frozenset("ACDEFGHIKLMNPQRSTVWY")
PROTEIN_AMBIGUITY = frozenset("BJOUXZ")
PROTEIN_ALPHABET = PROTEIN_LETTERS | PROTEIN_AMBIGUITY

#: Below this many letters, alphabet detection returns "unknown" instead of
#: guessing. Every single-letter amino-acid code is also a possible
#: nucleotide ambiguity code, so short strings are genuinely undecidable.
MIN_DETECTION_LENGTH = 4

_WHITESPACE = re.compile(r"\s+")

# Characters that routinely appear inside pasted sequences and carry no
# residue meaning: position numbers from GenBank/EMBL listings, alignment
# gaps, stop symbols and separators. They are silently ignored rather than
# reported as invalid, because rejecting a sequence copied out of a database
# record would be unhelpful and would not catch any real mistake.
IGNORABLE = frozenset(" \t\r\n0123456789-.*|/\\")


class SequenceError(ValueError):
    """Raised when input cannot be interpreted as the requested alphabet."""


@dataclass(frozen=True)
class FastaRecord:
    """One FASTA record. ``identifier`` is the first whitespace-delimited token."""

    identifier: str
    description: str
    sequence: str

    @property
    def header(self) -> str:
        return f"{self.identifier} {self.description}".strip()

    def __len__(self) -> int:  # pragma: no cover - trivial
        return len(self.sequence)


def strip_headers_and_whitespace(text: str) -> str:
    """Drop FASTA headers, digits, punctuation and whitespace; keep letters."""
    if not text:
        return ""
    kept: list[str] = []
    for line in str(text).splitlines():
        stripped = line.lstrip()
        if stripped.startswith((">", ";")):
            continue
        kept.append(line)
    return "".join(ch for ch in "".join(kept) if ch.isalpha())


def clean_dna(text: str, *, keep_ambiguous: bool = True) -> str:
    """Return an uppercase DNA string.

    ``U`` is folded to ``T`` so RNA input is accepted transparently. Symbols
    outside the IUPAC nucleotide alphabet are discarded; callers that need to
    reject them should use :func:`validate_dna`.
    """
    raw = strip_headers_and_whitespace(text).upper().replace("U", "T")
    allowed = DNA_ALPHABET if keep_ambiguous else DNA_BASES
    return "".join(ch for ch in raw if ch in allowed)


def clean_protein(text: str, *, keep_ambiguous: bool = True) -> str:
    """Return an uppercase amino-acid string."""
    raw = strip_headers_and_whitespace(text).upper()
    allowed = PROTEIN_ALPHABET if keep_ambiguous else PROTEIN_LETTERS
    return "".join(ch for ch in raw if ch in allowed)


def canonical_dna(text: str) -> str:
    """DNA restricted to A/C/G/T, used by calculations that need real bases."""
    return "".join(ch for ch in clean_dna(text) if ch in DNA_BASES)


def invalid_characters(text: str, alphabet: Iterable[str]) -> list[str]:
    """Characters in ``text`` that are neither in ``alphabet`` nor ignorable.

    The scan runs on the header-stripped text *before* non-letters are
    discarded, so a stray symbol such as ``@`` is reported rather than
    silently deleted. Digits, gaps and whitespace are ignored by design; see
    :data:`IGNORABLE`.
    """
    allowed = set(alphabet)
    seen: list[str] = []
    for line in str(text or "").splitlines():
        if line.lstrip().startswith((">", ";")):
            continue
        for ch in line.upper():
            if ch in allowed or ch in IGNORABLE:
                continue
            if ch not in seen:
                seen.append(ch)
    return seen


def validate_dna(text: str) -> str:
    """Strict DNA validation: raise :class:`SequenceError` on bad input."""
    if not str(text or "").strip():
        raise SequenceError("Empty input: paste a DNA sequence or load a FASTA file.")
    bad = invalid_characters(str(text).upper().replace("U", "T"), DNA_ALPHABET)
    if bad:
        raise SequenceError(
            "Input contains characters that are not IUPAC nucleotide codes: "
            + ", ".join(repr(b) for b in bad[:8])
        )
    seq = clean_dna(text)
    if not seq:
        raise SequenceError("No nucleotide letters found in the input.")
    return seq


def validate_protein(text: str) -> str:
    """Strict protein validation: raise :class:`SequenceError` on bad input."""
    if not str(text or "").strip():
        raise SequenceError("Empty input: paste a protein sequence or load a FASTA file.")
    bad = invalid_characters(text, PROTEIN_ALPHABET)
    if bad:
        raise SequenceError(
            "Input contains characters that are not amino-acid codes: "
            + ", ".join(repr(b) for b in bad[:8])
        )
    seq = clean_protein(text)
    if not seq:
        raise SequenceError("No amino-acid letters found in the input.")
    return seq


def looks_like_dna(text: str, threshold: float = 0.90) -> bool:
    """Heuristic: is this nucleotide text?

    A short peptide of only A/C/G/T-named residues is genuinely ambiguous, so
    sequences under 12 letters are not classified as DNA unless they contain
    no residue that is protein-only.
    """
    seq = strip_headers_and_whitespace(text).upper().replace("U", "T")
    if len(seq) < MIN_DETECTION_LENGTH:
        return False
    frac = sum(ch in DNA_ALPHABET for ch in seq) / len(seq)
    if frac < threshold:
        return False
    protein_only = PROTEIN_LETTERS - DNA_ALPHABET
    return not any(ch in protein_only for ch in seq)


def looks_like_protein(text: str, threshold: float = 0.90) -> bool:
    """Heuristic: is this amino-acid text?

    Like :func:`looks_like_dna`, this refuses to classify very short input.
    Under four letters there is genuinely no evidence either way, and
    guessing would send the user to the wrong tab with a confident-looking
    error message.
    """
    seq = strip_headers_and_whitespace(text).upper()
    if len(seq) < MIN_DETECTION_LENGTH:
        return False
    frac = sum(ch in PROTEIN_ALPHABET for ch in seq) / len(seq)
    return frac >= threshold and not looks_like_dna(text)


def detect_alphabet(text: str) -> str:
    """Return ``"dna"``, ``"protein"`` or ``"unknown"``."""
    if looks_like_dna(text):
        return "dna"
    if looks_like_protein(text):
        return "protein"
    return "unknown"


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------


def parse_fasta(text: str) -> list[FastaRecord]:
    """Parse multi-record FASTA text.

    Plain sequence text with no ``>`` header is accepted and returned as a
    single record called ``sequence_1``, which is what users paste most often.
    """
    records: list[FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        if header is None and not chunks:
            return
        raw = "".join(chunks)
        seq = "".join(ch for ch in raw if ch.isalpha()).upper()
        if header is None:
            ident, desc = f"sequence_{len(records) + 1}", ""
        else:
            parts = header.split(None, 1)
            ident = parts[0] if parts else f"sequence_{len(records) + 1}"
            desc = parts[1] if len(parts) > 1 else ""
        if seq:
            records.append(FastaRecord(ident, desc, seq))

    for line in str(text or "").splitlines():
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            chunks = []
        elif line.startswith(";"):
            continue
        else:
            chunks.append(line.strip())
    flush()
    return records


def iter_fasta_file(path: str) -> Iterator[FastaRecord]:
    """Stream records from a FASTA file without loading it all at once."""
    header: str | None = None
    chunks: list[str] = []
    index = 0

    def build() -> FastaRecord | None:
        nonlocal index
        raw = "".join(chunks)
        seq = "".join(ch for ch in raw if ch.isalpha()).upper()
        if not seq:
            return None
        index += 1
        if header is None:
            return FastaRecord(f"sequence_{index}", "", seq)
        parts = header.split(None, 1)
        return FastaRecord(
            parts[0] if parts else f"sequence_{index}",
            parts[1] if len(parts) > 1 else "",
            seq,
        )

    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                record = build()
                if record is not None:
                    yield record
                header = line[1:].strip()
                chunks = []
            elif line.startswith(";"):
                continue
            else:
                chunks.append(line.strip())
    record = build()
    if record is not None:
        yield record


def write_fasta(records: Iterable[FastaRecord], line_width: int = 60) -> str:
    """Serialise records to FASTA text with wrapped sequence lines."""
    out: list[str] = []
    for rec in records:
        out.append(f">{rec.header}".rstrip())
        seq = rec.sequence
        if line_width and line_width > 0:
            for i in range(0, len(seq), line_width):
                out.append(seq[i : i + line_width])
        else:
            out.append(seq)
    return "\n".join(out) + ("\n" if out else "")


__all__ = [
    "DNA_ALPHABET",
    "DNA_BASES",
    "IGNORABLE",
    "MIN_DETECTION_LENGTH",
    "PROTEIN_ALPHABET",
    "PROTEIN_LETTERS",
    "FastaRecord",
    "SequenceError",
    "canonical_dna",
    "clean_dna",
    "clean_protein",
    "detect_alphabet",
    "invalid_characters",
    "iter_fasta_file",
    "looks_like_dna",
    "looks_like_protein",
    "parse_fasta",
    "strip_headers_and_whitespace",
    "validate_dna",
    "validate_protein",
    "write_fasta",
]
