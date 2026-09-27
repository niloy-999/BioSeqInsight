"""Structured result objects shared by every layer of BioSeqInsight.

Every analytical or retrieval operation returns one of these dataclasses
rather than a preformatted string. The GUI, the CLI, the exporters and the
benchmark harness all consume the *same* objects, which is what makes the
benchmark numbers in the manuscript reproducible from the command line
without a display server.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from enum import Enum
from typing import Any


def _asdict(value: Any) -> Any:
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {k: _asdict(v) for k, v in dataclasses.asdict(value).items()}
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, dict):
        return {k: _asdict(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_asdict(v) for v in value]
    return value


class DictMixin:
    """Adds a JSON-friendly ``to_dict`` to dataclass results."""

    def to_dict(self) -> dict[str, Any]:
        return _asdict(self)


# --------------------------------------------------------------------------
# Sequence-level results
# --------------------------------------------------------------------------


@dataclass
class MotifHit(DictMixin):
    motif: str
    start: int  # 0-based, on the strand searched
    end: int
    strand: str = "+"


@dataclass
class ORF(DictMixin):
    strand: str
    frame: int
    start: int  # 0-based on the scanned strand
    end: int
    aa_length: int
    protein: str
    unterminated: bool = False


@dataclass
class TmResult(DictMixin):
    tm_c: float
    method: str
    length: int
    wallace_c: float | None = None
    na_mM: float | None = None
    primer_nM: float | None = None
    note: str | None = None


@dataclass
class SequenceResult(DictMixin):
    """Complete deterministic analysis of one nucleotide sequence."""

    identifier: str
    length: int
    alphabet: str  # "dna" or "rna"
    gc_percent: float
    composition: dict[str, int]
    composition_percent: dict[str, float]
    ambiguous_count: int
    reverse_complement: str | None = None
    mrna: str | None = None
    peptides: list[str] = field(default_factory=list)
    orfs: list[ORF] = field(default_factory=list)
    motifs: list[MotifHit] = field(default_factory=list)
    tm: TmResult | None = None
    warnings: list[str] = field(default_factory=list)
    software_version: str = ""


@dataclass
class SecondaryStructureSketch(DictMixin):
    length: int
    visual: str
    helix_percent: float
    sheet_percent: float
    coil_percent: float
    method: str


@dataclass
class ProteinResult(DictMixin):
    """Complete deterministic analysis of one amino-acid sequence."""

    identifier: str
    length: int
    molecular_weight: float
    gravy: float
    aromaticity: float
    isoelectric_point: float
    aliphatic_index: float
    extinction_coefficient_reduced: int
    extinction_coefficient_cystines: int
    composition: dict[str, int]
    composition_percent: dict[str, float]
    unknown_residues: int
    secondary_structure: SecondaryStructureSketch | None = None
    warnings: list[str] = field(default_factory=list)
    software_version: str = ""


# --------------------------------------------------------------------------
# Structure-level results
# --------------------------------------------------------------------------


class MappingLevel(str, Enum):
    """Formal sequence-to-structure mapping taxonomy (M0-M4).

    The same taxonomy is used by the software, the CLI, the exporters and the
    benchmark harness, so a reported mapping level always means exactly one
    thing.

    ``M0_NO_RESULT``
        No structural record was returned by any configured provider.
    ``M1_UNVERIFIED``
        A structure was returned but its biological identity was not
        established. Two distinct situations land here, and the accompanying
        message says which: either no comparison was possible (no query
        sequence, or no readable chain in the coordinates), or a comparison
        was made and the structure is clearly *not* the requested molecule.
        Neither case is evidence for the queried protein.
    ``M2_SEQUENCE_MATCH``
        The returned chain aligns to the query above the "related" thresholds
        but below the high-confidence thresholds. Typically a homologue, a
        different-species orthologue, or a fragment. Useful, but it is not
        the queried molecule.
    ``M3_HIGH_CONFIDENCE``
        Identity and coverage are both at or above the high-confidence
        thresholds, but the accession reported for the structure differs from
        the queried accession (or no accession was queried).
    ``M4_EXACT``
        The accession requested by the user is present in the structure's
        cross-references *and* the aligned identity and coverage reach the
        exact-match thresholds.
    """

    M0_NO_RESULT = "M0"
    M1_UNVERIFIED = "M1"
    M2_SEQUENCE_MATCH = "M2"
    M3_HIGH_CONFIDENCE = "M3"
    M4_EXACT = "M4"

    @property
    def label(self) -> str:
        return {
            "M0": "No result",
            "M1": "Structure found, identity unverified",
            "M2": "Sequence match (related record)",
            "M3": "High-confidence sequence match",
            "M4": "Exact accession and sequence match",
        }[self.value]

    @property
    def rank(self) -> int:
        return int(self.value[1:])


class RetrievalStatus(str, Enum):
    OK = "ok"
    NOT_FOUND = "not_found"
    SERVICE_ERROR = "service_error"
    TIMEOUT = "timeout"
    INVALID_INPUT = "invalid_input"
    SKIPPED = "skipped"


@dataclass
class ValidationReport(DictMixin):
    """Result of comparing a retrieved structure against the query."""

    mapping_level: MappingLevel
    identity_percent: float | None = None
    coverage_percent: float | None = None
    aligned_length: int | None = None
    query_length: int | None = None
    subject_length: int | None = None
    subject_chain: str | None = None
    accession_match: bool | None = None
    returned_accessions: list[str] = field(default_factory=list)
    message: str = ""

    @property
    def exact_biological_identity(self) -> bool:
        return self.mapping_level is MappingLevel.M4_EXACT


@dataclass
class AttemptRecord(DictMixin):
    """One HTTP attempt, retained so fault tolerance can be measured."""

    attempt: int
    url: str
    status_code: int | None
    error: str | None
    elapsed_s: float
    retried: bool = False


@dataclass
class StructureResult(DictMixin):
    """Unified structure record returned by every provider.

    ESM Atlas, AlphaFold DB and RCSB PDB return very different payloads. They
    are normalised here so that the GUI, the exporters and the benchmark can
    treat them identically and so that providers can be compared on equal
    terms.
    """

    source: str  # "esmatlas" | "alphafold" | "rcsb" | "local"
    query: str
    query_type: str  # "sequence" | "uniprot" | "pdb_id" | "file"
    status: RetrievalStatus
    identifier: str | None = None  # PDB ID, AF accession, or model name
    structure_url: str | None = None
    download_path: str | None = None
    sequence: str | None = None  # chain sequence parsed from coordinates
    confidence: float | None = None  # mean pLDDT where the source provides it
    confidence_kind: str | None = None  # "plddt" | "bfactor" | None
    experimental: bool = False
    response_time_s: float | None = None
    from_cache: bool = False
    attempts: list[AttemptRecord] = field(default_factory=list)
    cross_references: list[str] = field(default_factory=list)
    validation: ValidationReport | None = None
    error: str | None = None
    retrieved_at: str | None = None
    software_version: str = ""

    @property
    def ok(self) -> bool:
        return self.status is RetrievalStatus.OK

    @property
    def mapping_level(self) -> MappingLevel:
        if self.validation is not None:
            return self.validation.mapping_level
        return MappingLevel.M0_NO_RESULT if not self.ok else MappingLevel.M1_UNVERIFIED

    @property
    def recovered(self) -> bool:
        """True when the request eventually succeeded after >=1 real failure.

        Counts attempts that actually failed (an error, or an HTTP status of
        400 or above), not just the number of HTTP calls made. That
        distinction matters: some providers legitimately make more than one
        logical request on a completely clean run — AlphaFold DB always
        queries the prediction API to discover the current model URL before
        downloading it, which is two successful calls with nothing to
        recover from. Counting attempts alone would flag every ordinary
        AlphaFold retrieval as "recovered", which is what happened before
        this was fixed: a live run showed 67 of 67 successful AlphaFold
        retrievals marked as recovered-after-failure, none of which had
        actually failed at all.
        """
        if not self.ok:
            return False
        return any(
            attempt.error is not None
            or (attempt.status_code is not None and attempt.status_code >= 400)
            for attempt in self.attempts
        )


@dataclass
class BatchRow(DictMixin):
    """One row of a batch run: local analysis plus optional structure work."""

    identifier: str
    input_type: str
    length: int | None = None
    sequence: str | None = None  # retained for structure retrieval, not exported
    sequence_result: SequenceResult | None = None
    protein_result: ProteinResult | None = None
    structures: list[StructureResult] = field(default_factory=list)
    best_mapping: str = MappingLevel.M0_NO_RESULT.value
    elapsed_s: float = 0.0
    error: str | None = None


__all__ = [
    "ORF",
    "AttemptRecord",
    "BatchRow",
    "DictMixin",
    "MappingLevel",
    "MotifHit",
    "ProteinResult",
    "RetrievalStatus",
    "SecondaryStructureSketch",
    "SequenceResult",
    "StructureResult",
    "TmResult",
    "ValidationReport",
]
