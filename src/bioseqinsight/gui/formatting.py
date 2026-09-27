"""Result formatting for the graphical interface.

These functions turn result dataclasses into the text the user reads. They
are kept free of any Tkinter import so that the exact wording shown in the
application, including the warnings that prevent a non-exact structural match
from being read as a confirmed one, can be tested headless in continuous
integration.
"""

from __future__ import annotations

from ..models.results import (
    ProteinResult,
    SequenceResult,
    StructureResult,
)
from ..structures.manager import StructureSearchOutcome

DIVIDER = "-" * 62


def _wrap(text: str, width: int, indent: int = 0) -> str:
    """Wrap a bare sequence string at ``width`` characters."""
    pad = " " * indent
    return "\n".join(pad + text[i : i + width] for i in range(0, len(text), width)) or pad


def format_sequence_result(result: SequenceResult) -> str:
    lines = [
        f"Sequence analysis: {result.identifier}",
        DIVIDER,
        f"Length          {result.length} nt",
        f"GC content      {result.gc_percent:.2f}%",
        "Composition     "
        + "  ".join(
            f"{base}={count} ({result.composition_percent[base]:.1f}%)"
            for base, count in result.composition.items()
        ),
    ]
    if result.tm is not None:
        lines.append(f"Melting temp.   {result.tm.tm_c:.2f} C  [{result.tm.method}]")
        if result.tm.method == "nearest_neighbor":
            lines.append(
                f"                {result.tm.primer_nM:.0f} nM strand, "
                f"{result.tm.na_mM:.0f} mM Na+ (SantaLucia 1998)"
            )
        if result.tm.note:
            lines.append(f"                note: {result.tm.note}")
    lines.append(f"ORFs            {len(result.orfs)} found")
    for index, orf in enumerate(result.orfs[:10], 1):
        flag = "  [runs to end, no stop codon]" if orf.unterminated else ""
        lines.append(
            f"  {index:>2}. strand {orf.strand} frame {orf.frame}  "
            f"nt {orf.start}-{orf.end}  {orf.aa_length} aa{flag}"
        )
        lines.append(f"      {_wrap(orf.protein, 54, indent=6)}")
    if len(result.orfs) > 10:
        lines.append(f"  ... and {len(result.orfs) - 10} more")

    if result.motifs:
        lines += ["", f"Motif hits      {len(result.motifs)}"]
        for hit in result.motifs[:20]:
            lines.append(f"  {hit.motif}  {hit.strand} strand  {hit.start}-{hit.end}")
        if len(result.motifs) > 20:
            lines.append(f"  ... and {len(result.motifs) - 20} more")

    if result.warnings:
        lines += ["", "Warnings"]
        lines += [f"  - {w}" for w in result.warnings]
    lines += ["", f"BioSeqInsight {result.software_version}  (all values computed locally)"]
    return "\n".join(lines)


def format_protein_result(result: ProteinResult) -> str:
    lines = [
        f"Protein analysis: {result.identifier}",
        DIVIDER,
        f"Length              {result.length} aa",
        f"Molecular weight    {result.molecular_weight:.2f} Da "
        f"({result.molecular_weight / 1000:.3f} kDa)",
        f"GRAVY               {result.gravy:.4f} "
        f"({'hydrophobic' if result.gravy > 0 else 'hydrophilic'} overall)",
        f"Isoelectric point   {result.isoelectric_point:.2f}",
        f"Aromaticity         {result.aromaticity:.4f}",
        f"Aliphatic index     {result.aliphatic_index:.2f}",
        f"Extinction 280 nm   {result.extinction_coefficient_reduced} M-1 cm-1 (Cys reduced)",
        f"                    {result.extinction_coefficient_cystines} M-1 cm-1 (Cys as cystines)",
        "",
        "Amino-acid composition",
    ]
    present = [(aa, n) for aa, n in result.composition.items() if n]
    for index in range(0, len(present), 5):
        chunk = present[index : index + 5]
        lines.append(
            "  "
            + "  ".join(
                f"{aa} {n:>3} ({result.composition_percent[aa]:4.1f}%)" for aa, n in chunk
            )
        )

    sketch = result.secondary_structure
    if sketch is not None:
        lines += [
            "",
            "Secondary-structure propensity sketch",
            f"  {sketch.method}",
            f"  helix {sketch.helix_percent:.1f}%   sheet {sketch.sheet_percent:.1f}%   "
            f"coil {sketch.coil_percent:.1f}%",
            "  NOT a structure prediction; do not report these percentages as such.",
            "",
            _wrap(sketch.visual, 60, indent=2),
        ]
    if result.warnings:
        lines += ["", "Warnings"] + [f"  - {w}" for w in result.warnings]
    lines += ["", f"BioSeqInsight {result.software_version}  (all values computed locally)"]
    return "\n".join(lines)


def format_six_frames(frames: dict[str, str]) -> str:
    lines = ["Six-frame translation", DIVIDER]
    for name, peptide in frames.items():
        lines.append(f"Frame {name}  ({len(peptide)} aa)")
        lines.append(_wrap(peptide, 60, indent=2))
        lines.append("")
    lines.append("'*' marks a stop codon.")
    return "\n".join(lines)


def format_hydropathy(values: list[float], window: int) -> str:
    if not values:
        return "Sequence is shorter than the hydropathy window."
    width = 46
    low, high = min(values), max(values)
    span = (high - low) or 1.0
    lines = [
        f"Kyte-Doolittle hydropathy profile (window = {window})",
        DIVIDER,
        f"range {low:.2f} to {high:.2f}; positive values are hydrophobic",
        "",
    ]
    step = max(1, len(values) // 60)
    for position in range(0, len(values), step):
        value = values[position]
        filled = int((value - low) / span * width)
        marker = "#" * max(filled, 1)
        lines.append(f"{position + 1:>5} {value:>6.2f} |{marker}")
    return "\n".join(lines)




def format_structure_outcome(outcome: StructureSearchOutcome) -> str:
    """Render a full retrieval outcome, including every resource tried."""
    query = outcome.query
    preview = query.raw if len(query.raw) <= 60 else query.raw[:57] + "..."
    lines = [
        "Structure retrieval",
        DIVIDER,
        f"Query            {preview}",
        f"Interpreted as   {query.query_type}",
    ]
    if query.accession:
        lines.append(f"Accession        {query.accession}")
    if query.sequence:
        lines.append(f"Query length     {len(query.sequence)} aa")
    lines += [f"Elapsed          {outcome.elapsed_s:.2f} s", ""]

    for result in outcome.results:
        lines.append(f"[{result.source.upper()}] {result.status.value}")
        lines.extend(_format_result_body(result))
        lines.append("")

    best = outcome.best
    if best is not None:
        lines += [
            DIVIDER,
            f"BEST RESULT      {best.source} {best.identifier or ''}".rstrip(),
            f"MATCH STATUS     {outcome.best_mapping.value} - {outcome.best_mapping.label}",
        ]
        validation = best.validation
        if validation and validation.identity_percent is not None:
            lines += [
                f"Sequence identity {validation.identity_percent:.1f}%",
                f"Coverage          {validation.coverage_percent:.1f}%",
                f"Exact identity    {'YES' if validation.exact_biological_identity else 'NO'}",
            ]
        if best.download_path:
            lines.append(f"Coordinates       {best.download_path}")
    else:
        lines += [DIVIDER, "BEST RESULT      none", "MATCH STATUS     M0 - No result"]

    if outcome.warnings:
        lines += ["", "WARNINGS"]
        lines += [f"  ! {w}" for w in outcome.warnings]
    return "\n".join(lines)


def _format_result_body(result: StructureResult) -> list[str]:
    lines: list[str] = []
    if result.identifier:
        lines.append(f"  identifier     {result.identifier}")
    if result.structure_url:
        lines.append(f"  source         {result.structure_url}")
    if result.confidence is not None:
        label = "mean pLDDT" if result.confidence_kind == "plddt" else "mean B-factor"
        lines.append(f"  {label:<14} {result.confidence:.2f}")
    if result.attempts:
        detail = ", ".join(
            f"#{a.attempt} {a.status_code or 'network'}" + (" ok" if not a.error else "")
            for a in result.attempts
        )
        lines.append(f"  attempts       {len(result.attempts)} ({detail})")
        if result.recovered:
            lines.append("  recovery       succeeded after an initial failure")
    if result.from_cache:
        lines.append("  cache          served from the local cache")
    if result.response_time_s is not None:
        lines.append(f"  response time  {result.response_time_s:.2f} s")
    validation = result.validation
    if validation is not None:
        lines.append(
            f"  mapping        {validation.mapping_level.value} {validation.mapping_level.label}"
        )
        if validation.identity_percent is not None:
            lines.append(
                f"  identity       {validation.identity_percent:.1f}%   "
                f"coverage {validation.coverage_percent:.1f}%"
                + (f"   chain {validation.subject_chain}" if validation.subject_chain else "")
            )
        if validation.returned_accessions:
            lines.append(f"  cross-refs     {', '.join(validation.returned_accessions[:5])}")
        if validation.message:
            lines.append(f"  verdict        {validation.message}")
    if result.error:
        lines.append(f"  error          {result.error}")
    return lines


__all__ = [
    "DIVIDER",
    "format_hydropathy",
    "format_protein_result",
    "format_sequence_result",
    "format_six_frames",
    "format_structure_outcome",
]
