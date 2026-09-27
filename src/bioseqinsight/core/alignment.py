"""Pairwise alignment used to verify that a retrieved structure really is the
protein the user asked for.

This is the computational core of the v2.0 result-validation system. When a
structural resource returns a coordinate file, BioSeqInsight extracts the
chain sequence and aligns it against the query. Identity and coverage from
that alignment, combined with accession cross-references, determine the
mapping level (M0-M4) reported to the user.

Two algorithms are provided:

``needleman_wunsch``
    Full global alignment with affine-free linear gap penalties. Exact, and
    fast enough for the protein sizes involved (a 400x400 comparison takes a
    few hundred milliseconds in pure Python).
``semi_global``
    Global on the query, free end-gaps on the subject. This is the right
    model when a PDB chain is a fragment of the full-length UniProt sequence,
    which is extremely common and which a strict global alignment would
    wrongly penalise as low identity.

Deliberately *not* implemented: heuristic database search. BLAST-style search
is the job of the external services BioSeqInsight queries; this module only
verifies a single returned candidate.
"""

from __future__ import annotations

from dataclasses import dataclass

MATCH = 1
MISMATCH = -1
GAP = -2
_MAX_CELLS = 4_000_000  # ~2000x2000; guards against pathological inputs


@dataclass(frozen=True)
class AlignmentResult:
    """Outcome of one pairwise alignment."""

    score: int
    identities: int
    aligned_columns: int
    aligned_query_residues: int
    query_length: int
    subject_length: int
    aligned_query: str
    aligned_subject: str
    mode: str

    @property
    def identity_percent(self) -> float:
        """Identical columns as a percentage of aligned columns.

        Columns where both sequences have a residue are counted; end gaps in
        semi-global mode are excluded, because including them would report a
        perfectly matching 76-residue fragment of a 300-residue protein as
        25% identical.
        """
        if not self.aligned_columns:
            return 0.0
        return 100.0 * self.identities / self.aligned_columns

    @property
    def coverage_percent(self) -> float:
        """Fraction of the *query* covered by aligned subject residues."""
        if not self.query_length:
            return 0.0
        return 100.0 * self.aligned_query_residues / self.query_length


class AlignmentTooLarge(ValueError):
    """Raised when an alignment would exceed the safety limit on cells."""


def _check_size(a: str, b: str) -> None:
    if (len(a) + 1) * (len(b) + 1) > _MAX_CELLS:
        raise AlignmentTooLarge(
            f"Alignment of {len(a)} x {len(b)} residues exceeds the "
            f"{_MAX_CELLS}-cell safety limit. Compare a chain or domain instead."
        )


def _traceback_stats(aligned_a: str, aligned_b: str, mode: str) -> tuple[int, int, int]:
    """Return (identities, aligned columns, aligned query residues)."""
    identities = 0
    columns = 0
    covered = 0
    for x, y in zip(aligned_a, aligned_b):
        if x != "-" and y != "-":
            columns += 1
            covered += 1
            if x == y:
                identities += 1
        elif (x != "-" and y == "-") or (x == "-" and y != "-"):
            if mode == "global":
                columns += 1
    return identities, columns, covered


def _align(query: str, subject: str, free_end_gaps: bool) -> AlignmentResult:
    q = (query or "").upper()
    s = (subject or "").upper()
    _check_size(q, s)
    n, m = len(q), len(s)
    if n == 0 or m == 0:
        return AlignmentResult(
            score=0,
            identities=0,
            aligned_columns=0,
            aligned_query_residues=0,
            query_length=n,
            subject_length=m,
            aligned_query="-" * m if n == 0 else q,
            aligned_subject="-" * n if m == 0 else s,
            mode="semi_global" if free_end_gaps else "global",
        )

    # Dynamic-programming matrix stored as a flat list of rows.
    matrix: list[list[int]] = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        matrix[i][0] = matrix[i - 1][0] + GAP
    for j in range(1, m + 1):
        matrix[0][j] = 0 if free_end_gaps else matrix[0][j - 1] + GAP

    for i in range(1, n + 1):
        qi = q[i - 1]
        row = matrix[i]
        prev = matrix[i - 1]
        for j in range(1, m + 1):
            diag = prev[j - 1] + (MATCH if qi == s[j - 1] else MISMATCH)
            up = prev[j] + GAP
            left = row[j - 1] + GAP
            row[j] = diag if diag >= up and diag >= left else (up if up >= left else left)

    # Start of traceback.
    if free_end_gaps:
        best_j = max(range(m + 1), key=lambda j: matrix[n][j])
        i, j = n, best_j
        score = matrix[n][best_j]
        tail_q, tail_s = "-" * (m - best_j), s[best_j:]
    else:
        i, j = n, m
        score = matrix[n][m]
        tail_q = tail_s = ""

    out_q: list[str] = []
    out_s: list[str] = []
    while i > 0 and j > 0:
        current = matrix[i][j]
        if current == matrix[i - 1][j - 1] + (MATCH if q[i - 1] == s[j - 1] else MISMATCH):
            out_q.append(q[i - 1])
            out_s.append(s[j - 1])
            i -= 1
            j -= 1
        elif current == matrix[i - 1][j] + GAP:
            out_q.append(q[i - 1])
            out_s.append("-")
            i -= 1
        else:
            out_q.append("-")
            out_s.append(s[j - 1])
            j -= 1
    while i > 0:
        out_q.append(q[i - 1])
        out_s.append("-")
        i -= 1
    while j > 0:
        out_q.append("-")
        out_s.append(s[j - 1])
        j -= 1

    aligned_q = "".join(reversed(out_q)) + tail_q
    aligned_s = "".join(reversed(out_s)) + tail_s
    mode = "semi_global" if free_end_gaps else "global"
    identities, columns, covered = _traceback_stats(aligned_q, aligned_s, mode)
    return AlignmentResult(
        score=score,
        identities=identities,
        aligned_columns=columns,
        aligned_query_residues=covered,
        query_length=n,
        subject_length=m,
        aligned_query=aligned_q,
        aligned_subject=aligned_s,
        mode=mode,
    )


def needleman_wunsch(query: str, subject: str) -> AlignmentResult:
    """Exact global alignment with linear gap penalties."""
    return _align(query, subject, free_end_gaps=False)


def semi_global(query: str, subject: str) -> AlignmentResult:
    """Global on the query, free end gaps on the subject."""
    return _align(query, subject, free_end_gaps=True)


def quick_identity(query: str, subject: str) -> float:
    """Identity without alignment, for equal-length sequences.

    Used as a fast path: identical or equal-length sequences are extremely
    common in structure validation and do not need the full matrix.
    """
    q, s = (query or "").upper(), (subject or "").upper()
    if not q or not s or len(q) != len(s):
        return 0.0
    return 100.0 * sum(a == b for a, b in zip(q, s)) / len(q)


def compare(query: str, subject: str, allow_fragment: bool = True) -> AlignmentResult:
    """Validation entry point with a fast path for trivial cases."""
    q, s = (query or "").upper(), (subject or "").upper()
    if q and q == s:
        return AlignmentResult(
            score=MATCH * len(q),
            identities=len(q),
            aligned_columns=len(q),
            aligned_query_residues=len(q),
            query_length=len(q),
            subject_length=len(s),
            aligned_query=q,
            aligned_subject=s,
            mode="identical",
        )
    return semi_global(q, s) if allow_fragment else needleman_wunsch(q, s)


__all__ = [
    "AlignmentResult",
    "AlignmentTooLarge",
    "compare",
    "needleman_wunsch",
    "quick_identity",
    "semi_global",
]
