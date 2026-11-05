"""Result validation: does the returned structure actually correspond to the
protein the user asked for?

This module answers a gap in BioSeqInsight v1.0, whose benchmark reported a
structural hit for 25 of 25 queries but noted that for 7 of them the UniProt
accession attached to the hit differed from the query. v1.0 had no step that
compared the retrieved structure with the query sequence, so it could not say
whether such a hit was nonetheless a structure of the queried molecule (for
example a multi-chain entry that contains it), the same sequence filed under
another accession, or a genuine mismatch. v2.0 makes that comparison
explicitly.

v2.0 refuses to conflate "a structure was returned" with "this is your
protein". Every retrieved structure is aligned against the query sequence and
classified on an explicit five-level scale:

======  ====================================================================
Level   Meaning
======  ====================================================================
M0      No structural record was returned.
M1      A structure was returned but could not be compared (no query
        sequence, or no usable chain sequence in the coordinates).
M2      Aligns above the related thresholds: homologue, orthologue or
        fragment. Usable, but it is not the queried molecule.
M3      Identity and coverage reach the high-confidence thresholds, but the
        accession was not confirmed.
M4      The queried accession appears in the structure's cross-references
        and identity and coverage reach the exact thresholds.
======  ====================================================================

Thresholds live in :class:`~bioseqinsight.config.settings.Settings` so that a
user can tighten or relax them and re-run the benchmark, and the values used
are recorded in every exported result.
"""

from __future__ import annotations

from ..config.settings import Settings
from ..core.alignment import AlignmentTooLarge, compare
from ..models.results import MappingLevel, ValidationReport
from .pdbio import ParsedStructure, find_accessions, parse_pdb


def classify(
    identity: float | None,
    coverage: float | None,
    accession_match: bool | None,
    settings: Settings,
) -> MappingLevel:
    """Map identity, coverage and accession agreement onto M0-M4."""
    if identity is None or coverage is None:
        return MappingLevel.M1_UNVERIFIED
    if (
        accession_match
        and identity >= settings.exact_identity_threshold
        and coverage >= settings.exact_coverage_threshold
    ):
        return MappingLevel.M4_EXACT
    if (
        identity >= settings.high_identity_threshold
        and coverage >= settings.high_coverage_threshold
    ):
        return MappingLevel.M3_HIGH_CONFIDENCE
    if (
        identity >= settings.related_identity_threshold
        and coverage >= settings.related_coverage_threshold
    ):
        return MappingLevel.M2_SEQUENCE_MATCH
    # The comparison succeeded but the structure falls below even the
    # "related" thresholds: it is not the requested molecule and must not be
    # promoted to M2, which means "related record".
    return MappingLevel.M1_UNVERIFIED


def is_partial_structure(
    identity: float | None,
    coverage: float | None,
    accession_match: bool | None,
    settings: Settings,
) -> bool:
    """True for a confirmed but incomplete structure of the queried protein.

    The accession is confirmed by the structure's cross-references and the
    aligned region is at or above the high-identity threshold, but the
    structure covers less of the query than the high-coverage threshold
    (a domain construct, for example). Such a result is a real structure of
    the requested protein, so telling the user it "does not correspond" to
    it would be wrong; it is also not full-length, so it is never M4.

    This is a qualifier only: it does not change the mapping level.
    """
    if identity is None or coverage is None or not accession_match:
        return False
    return (
        identity >= settings.high_identity_threshold
        and coverage < settings.high_coverage_threshold
    )


def validate_structure(
    pdb_text: str,
    *,
    query_sequence: str | None,
    query_accession: str | None,
    settings: Settings,
    extra_accessions: list[str] | None = None,
    parsed: ParsedStructure | None = None,
) -> ValidationReport:
    """Compare a coordinate file against the query and classify the match.

    ``extra_accessions`` carries cross-references obtained from the provider's
    metadata API (for example the RCSB polymer-entity record), which are often
    richer than the DBREF lines inside the coordinate file itself.
    """
    structure = parsed if parsed is not None else parse_pdb(pdb_text)
    returned = list(structure.accessions)
    for accession in extra_accessions or []:
        if accession and accession not in returned:
            returned.append(accession)

    accession_match: bool | None = None
    if query_accession:
        wanted = query_accession.strip().upper()
        accession_match = wanted in {a.upper() for a in returned} if returned else False

    if not query_sequence or not structure.chains:
        reason = (
            "no query sequence was supplied for comparison"
            if not query_sequence
            else "the coordinate file contained no readable polymer chain"
        )
        return ValidationReport(
            mapping_level=MappingLevel.M1_UNVERIFIED,
            accession_match=accession_match,
            returned_accessions=returned,
            subject_length=(
                structure.longest_chain.residue_count if structure.longest_chain else None
            ),
            query_length=len(query_sequence) if query_sequence else None,
            message=(
                "Structure retrieved, but biological identity is unverified because "
                f"{reason}. Do not assume this is the requested molecule."
            ),
        )

    chain = structure.best_matching_chain(query_sequence)
    assert chain is not None  # guaranteed by the emptiness check above
    try:
        alignment = compare(query_sequence, chain.sequence, allow_fragment=True)
    except AlignmentTooLarge as exc:
        return ValidationReport(
            mapping_level=MappingLevel.M1_UNVERIFIED,
            accession_match=accession_match,
            returned_accessions=returned,
            query_length=len(query_sequence),
            subject_length=chain.residue_count,
            subject_chain=chain.chain_id,
            message=f"Identity could not be verified: {exc}",
        )

    identity = alignment.identity_percent
    coverage = alignment.coverage_percent
    level = classify(identity, coverage, accession_match, settings)
    partial = is_partial_structure(identity, coverage, accession_match, settings)
    return ValidationReport(
        mapping_level=level,
        identity_percent=identity,
        coverage_percent=coverage,
        aligned_length=alignment.aligned_columns,
        query_length=alignment.query_length,
        subject_length=alignment.subject_length,
        subject_chain=chain.chain_id,
        accession_match=accession_match,
        returned_accessions=returned,
        partial_structure=partial,
        message=describe(
            level, identity, coverage, accession_match, query_accession, returned, partial=partial
        ),
    )


def describe(
    level: MappingLevel,
    identity: float | None,
    coverage: float | None,
    accession_match: bool | None,
    query_accession: str | None,
    returned: list[str],
    partial: bool = False,
) -> str:
    """Human-readable verdict shown in the GUI and written to exports."""
    if level is MappingLevel.M0_NO_RESULT:
        return "No structural record was returned by any configured resource."
    if partial and level in (MappingLevel.M1_UNVERIFIED, MappingLevel.M2_SEQUENCE_MATCH):
        return (
            f"Partial structure of the requested protein: the structure cross-references "
            f"{query_accession} and matches the query at {identity:.1f}% identity over the "
            f"region it contains, but it covers only {coverage:.1f}% of your sequence. "
            "It is a real structure of this protein, not a full-length one."
        )
    if level is MappingLevel.M1_UNVERIFIED:
        if identity is None:
            return "Structure retrieved; biological identity could not be verified."
        return (
            f"WARNING: the retrieved structure does not correspond to the requested "
            f"protein (identity {identity:.1f}%, coverage {coverage:.1f}%), which is "
            "below the threshold for even a related record."
            + (f" It is associated with {', '.join(returned[:3])}." if returned else "")
        )

    numbers = f"identity {identity:.1f}%, coverage {coverage:.1f}%"
    if level is MappingLevel.M4_EXACT:
        return (
            f"Exact biological identity confirmed ({numbers}); the structure "
            f"cross-references {query_accession}."
        )
    if level is MappingLevel.M3_HIGH_CONFIDENCE:
        detail = (
            f"but the structure reports {', '.join(returned[:3])} rather than {query_accession}"
            if query_accession and returned and accession_match is False
            else "and no conflicting accession was reported"
        )
        return (
            f"High-confidence sequence match ({numbers}), {detail}. "
            "Treat this as the same or a closely related molecule, not as a confirmed identity."
        )
    return (
        f"WARNING: the returned structure may not correspond to the requested protein "
        f"({numbers})."
        + (
            f" The structure is associated with {', '.join(returned[:3])}."
            if returned
            else ""
        )
    )


def accessions_from_text(text: str) -> list[str]:
    """Expose accession scraping for provider metadata payloads."""
    return find_accessions(text)


def no_result_report(message: str = "") -> ValidationReport:
    return ValidationReport(
        mapping_level=MappingLevel.M0_NO_RESULT,
        message=message or "No structural record was returned by any configured resource.",
    )


__all__ = [
    "accessions_from_text",
    "classify",
    "describe",
    "is_partial_structure",
    "no_result_report",
    "validate_structure",
]
