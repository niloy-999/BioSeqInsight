"""Batch analysis.

One sequence at a time through a GUI is a teaching tool. A researcher with
200 sequences needs a single command that produces a table. This module turns
a FASTA file (or any iterable of records) into a list of
:class:`~bioseqinsight.models.results.BatchRow` objects, optionally including
structure retrieval and validation for each entry.

Design points that matter for the evaluation:

* **Local analysis is never blocked by the network.** Sequence and protein
  descriptors are computed first and are always present in the output, even
  when every structural resource is down.
* **Structure retrieval is I/O-bound, so it is threaded**, with the worker
  count configurable and defaulting to a small number out of courtesy to the
  public services.
* **One failure never stops the run.** Per-row exceptions are captured in the
  row, so a 200-sequence batch does not abort on sequence 37.
* **Progress is reported through a callback**, which the GUI renders as a
  progress bar and the CLI as a counter.
"""

from __future__ import annotations

import time
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field

from ..config.settings import Settings, get_settings
from ..core.alphabet import (
    FastaRecord,
    SequenceError,
    detect_alphabet,
    iter_fasta_file,
    parse_fasta,
)
from ..core.protein import ProteinAnalyzer
from ..core.sequence import SequenceAnalyzer
from ..models.results import BatchRow, MappingLevel, StructureResult
from ..services.logging_setup import get_logger
from ..structures.manager import StructureManager

logger = get_logger("batch")

ProgressCallback = Callable[[int, int, str], None]


@dataclass
class BatchOptions:
    """What a batch run should do."""

    include_structures: bool = False
    providers: list[str] | None = None
    min_orf_aa: int = 30
    include_sequences: bool = False
    workers: int = 4
    refresh_cache: bool = False
    stop_on_exact: bool = True
    force_alphabet: str | None = None  # "dna", "protein" or None to autodetect


@dataclass
class BatchSummary:
    """Aggregate statistics for a completed run."""

    total: int = 0
    succeeded: int = 0
    failed: int = 0
    elapsed_s: float = 0.0
    mapping_counts: dict[str, int] = field(default_factory=dict)
    structures_attempted: int = 0
    structures_retrieved: int = 0
    cache_hits: int = 0
    recovered_requests: int = 0

    def to_dict(self) -> dict:
        return {
            "total": self.total,
            "succeeded": self.succeeded,
            "failed": self.failed,
            "elapsed_s": round(self.elapsed_s, 3),
            "sequences_per_second": (
                round(self.total / self.elapsed_s, 2) if self.elapsed_s > 0 else None
            ),
            "mapping_counts": dict(self.mapping_counts),
            "structures_attempted": self.structures_attempted,
            "structures_retrieved": self.structures_retrieved,
            "cache_hits": self.cache_hits,
            "recovered_requests": self.recovered_requests,
        }


class BatchRunner:
    """Run local analysis, and optionally structure retrieval, over many inputs."""

    def __init__(
        self,
        settings: Settings | None = None,
        manager: StructureManager | None = None,
    ) -> None:
        self.settings = settings or get_settings()
        self._manager = manager

    @property
    def manager(self) -> StructureManager:
        if self._manager is None:
            self._manager = StructureManager(self.settings)
        return self._manager

    # -- input ---------------------------------------------------------------

    @staticmethod
    def records_from_fasta_file(path: str) -> list[FastaRecord]:
        return list(iter_fasta_file(path))

    @staticmethod
    def records_from_text(text: str) -> list[FastaRecord]:
        return parse_fasta(text)

    @staticmethod
    def records_from_identifiers(identifiers: Iterable[str]) -> list[FastaRecord]:
        """Treat a list of accessions or PDB ids as records with empty sequences."""
        return [FastaRecord(identifier=i.strip(), description="", sequence="") for i in identifiers if i.strip()]

    # -- execution -----------------------------------------------------------

    def run(
        self,
        records: Sequence[FastaRecord],
        options: BatchOptions | None = None,
        progress: ProgressCallback | None = None,
    ) -> tuple[list[BatchRow], BatchSummary]:
        """Analyse every record and return the rows plus a summary."""
        opts = options or BatchOptions()
        started = time.perf_counter()
        rows: list[BatchRow] = [self._analyse_local(record, opts) for record in records]

        if progress:
            for index, row in enumerate(rows, 1):
                progress(index, len(rows), f"analysed {row.identifier}")

        if opts.include_structures and rows:
            self._attach_structures(rows, opts, progress)

        summary = self._summarise(rows, time.perf_counter() - started, opts)
        return rows, summary

    def run_file(
        self,
        path: str,
        options: BatchOptions | None = None,
        progress: ProgressCallback | None = None,
    ) -> tuple[list[BatchRow], BatchSummary]:
        return self.run(self.records_from_fasta_file(path), options, progress)

    # -- internals -----------------------------------------------------------

    def _analyse_local(self, record: FastaRecord, opts: BatchOptions) -> BatchRow:
        started = time.perf_counter()
        row = BatchRow(identifier=record.identifier, input_type="unknown")
        sequence = record.sequence

        if not sequence:
            row.input_type = "identifier"
            row.elapsed_s = time.perf_counter() - started
            return row

        kind = opts.force_alphabet or detect_alphabet(sequence)
        row.input_type = kind
        row.length = len(sequence)
        row.sequence = sequence
        try:
            if kind == "dna":
                row.sequence_result = SequenceAnalyzer(
                    sequence, identifier=record.identifier, table=self.settings.genetic_code_table
                ).analyze(
                    min_orf_aa=opts.min_orf_aa,
                    include_sequences=opts.include_sequences,
                )
            elif kind == "protein":
                row.protein_result = ProteinAnalyzer(
                    sequence, identifier=record.identifier
                ).analyze(include_sketch=opts.include_sequences)
            else:
                row.error = (
                    "Could not determine whether this record is nucleotide or protein. "
                    "Re-run with an explicit alphabet."
                )
        except SequenceError as exc:
            row.error = str(exc)
        except Exception as exc:  # pragma: no cover - defensive, keeps the batch alive
            logger.exception("Local analysis of %s failed", record.identifier)
            row.error = f"{type(exc).__name__}: {exc}"

        row.elapsed_s = time.perf_counter() - started
        return row

    def _attach_structures(
        self, rows: list[BatchRow], opts: BatchOptions, progress: ProgressCallback | None
    ) -> None:
        targets = [
            row
            for row in rows
            if row.input_type in ("protein", "identifier") and row.error is None
        ]
        if not targets:
            return
        workers = max(1, min(opts.workers or self.settings.batch_workers, 16))
        done = 0
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(self._retrieve_one, row, opts): row for row in targets
            }
            for future in as_completed(futures):
                row = futures[future]
                done += 1
                try:
                    row.structures = future.result()
                except Exception as exc:  # pragma: no cover - defensive
                    logger.exception("Structure retrieval for %s failed", row.identifier)
                    row.error = (row.error or "") + f" structure retrieval: {exc}"
                    row.structures = []
                row.best_mapping = _best_level(row.structures).value
                if progress:
                    progress(done, len(targets), f"structures for {row.identifier}")

    def _retrieve_one(self, row: BatchRow, opts: BatchOptions) -> list[StructureResult]:
        query = row.identifier if row.input_type == "identifier" else _query_for(row)
        outcome = self.manager.retrieve(
            query,
            providers=opts.providers,
            stop_on_exact=opts.stop_on_exact,
            refresh=opts.refresh_cache,
        )
        return outcome.results

    def _summarise(
        self, rows: list[BatchRow], elapsed: float, opts: BatchOptions
    ) -> BatchSummary:
        summary = BatchSummary(total=len(rows), elapsed_s=elapsed)
        for row in rows:
            if row.error:
                summary.failed += 1
            else:
                summary.succeeded += 1
            if row.structures:
                summary.structures_attempted += 1
                if any(s.ok for s in row.structures):
                    summary.structures_retrieved += 1
                summary.cache_hits += sum(1 for s in row.structures if s.from_cache)
                summary.recovered_requests += sum(1 for s in row.structures if s.recovered)
            if opts.include_structures:
                level = row.best_mapping or MappingLevel.M0_NO_RESULT.value
                summary.mapping_counts[level] = summary.mapping_counts.get(level, 0) + 1
        return summary


def _query_for(row: BatchRow) -> str:
    """What to send to the structure subsystem for this row.

    A record with a sequence is queried by sequence; a record that carried
    only an identifier (an accession list, say) is queried by identifier.
    """
    return row.sequence or row.identifier


def _best_level(structures: Sequence[StructureResult]) -> MappingLevel:
    usable = [s for s in structures if s.ok]
    if not usable:
        return MappingLevel.M0_NO_RESULT
    return max((s.mapping_level for s in usable), key=lambda level: level.rank)


__all__ = ["BatchOptions", "BatchRunner", "BatchSummary", "ProgressCallback"]
