"""Exporters.

Everything BioSeqInsight computes can leave the application in a
machine-readable form. A researcher can therefore analyse 200 sequences and
hand the resulting table to R, pandas or a spreadsheet without retyping
anything, which is the practical difference between a teaching demonstration
and research software.

Formats: CSV and TSV for tables, JSON for full nested results (including
attempt histories and validation reports), FASTA for sequences, and a
self-contained HTML report for sharing.
"""

from __future__ import annotations

import csv
import html
import io
import json
import os
from collections.abc import Iterable, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .. import __version__
from ..core.alphabet import FastaRecord, write_fasta
from ..models.results import BatchRow, MappingLevel, ProteinResult, SequenceResult

BATCH_COLUMNS: tuple[str, ...] = (
    "identifier",
    "input_type",
    "length",
    "gc_percent",
    "tm_c",
    "tm_method",
    "orf_count",
    "longest_orf_aa",
    "molecular_weight",
    "gravy",
    "isoelectric_point",
    "aromaticity",
    "aliphatic_index",
    "structure_source",
    "structure_id",
    "structure_status",
    "mapping_level",
    "mapping_label",
    "identity_percent",
    "coverage_percent",
    "confidence",
    "confidence_kind",
    "from_cache",
    "attempts",
    "response_time_s",
    "elapsed_s",
    "error",
)


def _round(value: Any, digits: int = 4) -> Any:
    return round(value, digits) if isinstance(value, float) else value


def batch_row_to_record(row: BatchRow) -> dict[str, Any]:
    """Flatten one :class:`BatchRow` into a single table row."""
    record: dict[str, Any] = dict.fromkeys(BATCH_COLUMNS, "")
    record["identifier"] = row.identifier
    record["input_type"] = row.input_type
    record["length"] = row.length if row.length is not None else ""
    record["elapsed_s"] = _round(row.elapsed_s)
    record["error"] = row.error or ""

    if row.sequence_result is not None:
        seq = row.sequence_result
        record["gc_percent"] = _round(seq.gc_percent, 2)
        record["orf_count"] = len(seq.orfs)
        record["longest_orf_aa"] = seq.orfs[0].aa_length if seq.orfs else 0
        if seq.tm is not None:
            record["tm_c"] = _round(seq.tm.tm_c, 2)
            record["tm_method"] = seq.tm.method

    if row.protein_result is not None:
        prot = row.protein_result
        record["molecular_weight"] = _round(prot.molecular_weight, 2)
        record["gravy"] = _round(prot.gravy, 4)
        record["isoelectric_point"] = _round(prot.isoelectric_point, 2)
        record["aromaticity"] = _round(prot.aromaticity, 4)
        record["aliphatic_index"] = _round(prot.aliphatic_index, 2)

    best = _best_structure(row)
    if best is not None:
        record["structure_source"] = best.source
        record["structure_id"] = best.identifier or ""
        record["structure_status"] = best.status.value
        record["mapping_level"] = best.mapping_level.value
        record["mapping_label"] = best.mapping_level.label
        record["confidence"] = _round(best.confidence, 2)
        record["confidence_kind"] = best.confidence_kind or ""
        record["from_cache"] = best.from_cache
        record["attempts"] = len(best.attempts)
        record["response_time_s"] = _round(best.response_time_s)
        if best.validation is not None:
            record["identity_percent"] = _round(best.validation.identity_percent, 2)
            record["coverage_percent"] = _round(best.validation.coverage_percent, 2)
    else:
        record["mapping_level"] = MappingLevel.M0_NO_RESULT.value
        record["mapping_label"] = MappingLevel.M0_NO_RESULT.label
    return record


def _best_structure(row: BatchRow):
    usable = [s for s in row.structures if s.ok]
    if not usable:
        return row.structures[0] if row.structures else None
    return max(usable, key=lambda s: s.mapping_level.rank)


def rows_to_delimited(rows: Sequence[BatchRow], delimiter: str = ",") -> str:
    """Serialise batch rows as delimited text."""
    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer, fieldnames=list(BATCH_COLUMNS), delimiter=delimiter, lineterminator="\n"
    )
    writer.writeheader()
    for row in rows:
        writer.writerow(batch_row_to_record(row))
    return buffer.getvalue()


def write_csv(rows: Sequence[BatchRow], path: str | os.PathLike[str]) -> str:
    return _write(path, rows_to_delimited(rows, ","))


def write_tsv(rows: Sequence[BatchRow], path: str | os.PathLike[str]) -> str:
    return _write(path, rows_to_delimited(rows, "\t"))


def write_json(payload: Any, path: str | os.PathLike[str], indent: int = 2) -> str:
    return _write(path, to_json(payload, indent=indent))


def to_json(payload: Any, indent: int = 2) -> str:
    """Serialise results (or anything containing them) to JSON."""
    return json.dumps(_jsonable(payload), indent=indent, sort_keys=False, default=str)


def _jsonable(value: Any) -> Any:
    if hasattr(value, "to_dict"):
        return value.to_dict()
    if isinstance(value, dict):
        return {k: _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def write_sequences_fasta(
    records: Iterable[FastaRecord], path: str | os.PathLike[str]
) -> str:
    return _write(path, write_fasta(records))


def _write(path: str | os.PathLike[str], text: str) -> str:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text, encoding="utf-8")
    return str(target)


# --------------------------------------------------------------------------
# HTML report
# --------------------------------------------------------------------------

_REPORT_CSS = """
:root { color-scheme: light dark; --fg:#1f2328; --bg:#ffffff; --muted:#656d76;
  --line:#d0d7de; --accent:#0969da; --warn:#9a6700; --bad:#cf222e; --good:#1a7f37; }
@media (prefers-color-scheme: dark) {
  :root { --fg:#e6edf3; --bg:#0d1117; --muted:#8b949e; --line:#30363d;
    --accent:#4493f8; --warn:#d29922; --bad:#f85149; --good:#3fb950; } }
* { box-sizing: border-box; }
body { margin:0; padding:32px; max-width:1100px; margin-inline:auto; color:var(--fg);
  background:var(--bg); font-family: ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif;
  line-height:1.55; }
h1 { font-size:1.6rem; margin:0 0 4px; } h2 { font-size:1.15rem; margin:32px 0 8px;
  border-bottom:1px solid var(--line); padding-bottom:6px; }
.meta { color:var(--muted); font-size:.85rem; margin-bottom:24px; }
table { border-collapse:collapse; width:100%; font-size:.85rem; }
th, td { text-align:left; padding:6px 10px; border-bottom:1px solid var(--line);
  white-space:nowrap; }
th { font-weight:600; color:var(--muted); text-transform:uppercase; font-size:.7rem;
  letter-spacing:.04em; }
.wrap { overflow-x:auto; }
.M4 { color:var(--good); font-weight:600; } .M3 { color:var(--accent); }
.M2 { color:var(--warn); } .M1, .M0 { color:var(--bad); }
.cards { display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:12px; }
.card { border:1px solid var(--line); border-radius:8px; padding:12px 14px; }
.card .n { font-size:1.4rem; font-weight:650; } .card .l { color:var(--muted); font-size:.75rem; }
footer { margin-top:40px; color:var(--muted); font-size:.75rem;
  border-top:1px solid var(--line); padding-top:12px; }
"""


def build_html_report(
    rows: Sequence[BatchRow],
    title: str = "BioSeqInsight analysis report",
    notes: str = "",
) -> str:
    """Build a self-contained HTML report for a batch run."""
    records = [batch_row_to_record(row) for row in rows]
    counts: dict[str, int] = {}
    for record in records:
        counts[record["mapping_level"] or "M0"] = counts.get(record["mapping_level"] or "M0", 0) + 1

    cards = "".join(
        f'<div class="card"><div class="n">{value}</div><div class="l">{html.escape(label)}</div></div>'
        for label, value in (
            ("Sequences analysed", len(records)),
            ("Exact matches (M4)", counts.get("M4", 0)),
            ("High confidence (M3)", counts.get("M3", 0)),
            ("Related only (M2)", counts.get("M2", 0)),
            ("Unverified (M1)", counts.get("M1", 0)),
            ("No structure (M0)", counts.get("M0", 0)),
        )
    )

    header = "".join(f"<th>{html.escape(c.replace('_', ' '))}</th>" for c in BATCH_COLUMNS)
    body_rows = []
    for record in records:
        cells = []
        for column in BATCH_COLUMNS:
            value = html.escape(str(record[column]))
            css = f' class="{record["mapping_level"]}"' if column == "mapping_level" else ""
            cells.append(f"<td{css}>{value}</td>")
        body_rows.append("<tr>" + "".join(cells) + "</tr>")

    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    note_block = f"<p>{html.escape(notes)}</p>" if notes else ""
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title><style>{_REPORT_CSS}</style></head>
<body>
<h1>{html.escape(title)}</h1>
<div class="meta">Generated {generated} by BioSeqInsight {__version__}</div>
{note_block}
<h2>Summary</h2>
<div class="cards">{cards}</div>
<h2>Mapping levels</h2>
<p style="font-size:.85rem">
M4 exact accession and sequence match &middot; M3 high-confidence sequence match &middot;
M2 related record only &middot; M1 structure retrieved but unverified &middot; M0 no structure.
A result below M4 is not evidence that the retrieved coordinates represent the queried molecule.
</p>
<h2>Results</h2>
<div class="wrap"><table><thead><tr>{header}</tr></thead>
<tbody>{"".join(body_rows)}</tbody></table></div>
<footer>
Structures were retrieved from third-party public resources (RCSB PDB, AlphaFold DB, ESM Atlas).
BioSeqInsight does not generate structural predictions itself; cite the underlying resources.
</footer>
</body></html>
"""


def write_html_report(
    rows: Sequence[BatchRow], path: str | os.PathLike[str], **kwargs
) -> str:
    return _write(path, build_html_report(rows, **kwargs))


def sequence_result_to_record(result: SequenceResult) -> dict[str, Any]:
    """Flatten a single sequence result for one-row export."""
    return {
        "identifier": result.identifier,
        "length": result.length,
        "gc_percent": _round(result.gc_percent, 2),
        **{f"count_{base}": count for base, count in result.composition.items()},
        "tm_c": _round(result.tm.tm_c, 2) if result.tm else "",
        "tm_method": result.tm.method if result.tm else "",
        "orf_count": len(result.orfs),
        "longest_orf_aa": result.orfs[0].aa_length if result.orfs else 0,
        "software_version": result.software_version,
    }


def protein_result_to_record(result: ProteinResult) -> dict[str, Any]:
    """Flatten a single protein result for one-row export."""
    return {
        "identifier": result.identifier,
        "length": result.length,
        "molecular_weight": _round(result.molecular_weight, 2),
        "gravy": _round(result.gravy, 4),
        "isoelectric_point": _round(result.isoelectric_point, 2),
        "aromaticity": _round(result.aromaticity, 4),
        "aliphatic_index": _round(result.aliphatic_index, 2),
        "extinction_reduced": result.extinction_coefficient_reduced,
        "extinction_cystines": result.extinction_coefficient_cystines,
        **{f"aa_{aa}": count for aa, count in result.composition.items()},
        "software_version": result.software_version,
    }


def write_records_csv(
    records: Sequence[dict[str, Any]], path: str | os.PathLike[str], delimiter: str = ","
) -> str:
    """Write arbitrary flat records, using the union of their keys as columns."""
    if not records:
        return _write(path, "")
    columns: list[str] = []
    for record in records:
        for key in record:
            if key not in columns:
                columns.append(key)
    buffer = io.StringIO()
    writer = csv.DictWriter(
        buffer, fieldnames=columns, delimiter=delimiter, lineterminator="\n", restval=""
    )
    writer.writeheader()
    writer.writerows(records)
    return _write(path, buffer.getvalue())


__all__ = [
    "BATCH_COLUMNS",
    "batch_row_to_record",
    "build_html_report",
    "protein_result_to_record",
    "rows_to_delimited",
    "sequence_result_to_record",
    "to_json",
    "write_csv",
    "write_html_report",
    "write_json",
    "write_records_csv",
    "write_sequences_fasta",
    "write_tsv",
]
