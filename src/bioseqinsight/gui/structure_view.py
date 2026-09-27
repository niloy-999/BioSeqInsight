"""Structure-retrieval and batch tabs.

The structure tab is where v2.0 differs most visibly from v1.0. Instead of
printing "structure downloaded", it shows every resource that was queried,
how many attempts each took, and an explicit mapping level with the identity
and coverage behind it. A result that is not an exact match says so in the
interface, not only in the manuscript.
"""

from __future__ import annotations

import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

from ..core.alphabet import parse_fasta
from ..io.export import to_json, write_csv, write_html_report, write_json, write_tsv
from ..structures.manager import StructureSearchOutcome
from ..structures.visualization import (
    open_in_browser,
    open_in_external_viewer,
    write_viewer_html,
)
from ..workflows.batch import BatchOptions
from .formatting import DIVIDER, format_structure_outcome
from .widgets import BackgroundRunner, LabelledText, MappingBadge, bind_copy


class StructureView(ttk.Frame):
    """The Structure Retrieval tab."""

    def __init__(self, master, app, **kwargs):
        super().__init__(master, padding=12, **kwargs)
        self.app = app
        self.runner = BackgroundRunner(self)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)

        self.input = LabelledText(
            self,
            "Protein sequence, UniProt accession (P24941) or PDB identifier (1AQ1)",
            height=6,
        )
        self.input.grid(row=0, column=0, columnspan=2, sticky="nsew", pady=(0, 10))

        controls = ttk.Frame(self)
        controls.grid(row=1, column=0, sticky="nw", padx=(0, 12))
        right = ttk.Frame(self)
        right.grid(row=1, column=1, sticky="nsew")
        right.columnconfigure(0, weight=1)
        right.rowconfigure(1, weight=1)

        self.badge = MappingBadge(right)
        self.badge.grid(row=0, column=0, sticky="w", pady=(0, 4))
        self.output = LabelledText(right, "Result", height=18, readonly=True)
        self.output.grid(row=1, column=0, sticky="nsew")
        bind_copy(self.output.text)

        actions = [
            ("Open FASTA...", self.open_fasta),
            ("Retrieve and validate", self.run_retrieve),
            ("Query all resources", self.run_retrieve_all),
            ("AlphaFold DB only", lambda: self.run_retrieve(providers=["alphafold"])),
            ("RCSB PDB only", lambda: self.run_retrieve(providers=["rcsb"])),
            ("ESM Atlas fold only", lambda: self.run_retrieve(providers=["esmatlas"])),
            ("Load local PDB...", self.load_local),
            ("View in 3Dmol.js", self.view_browser),
            ("Open external viewer", self.view_external),
        ]
        for index, (label, command) in enumerate(actions):
            ttk.Button(controls, text=label, width=22, command=command).grid(
                row=index, column=0, sticky="ew", pady=2
            )

        self.refresh_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            controls, text="Ignore cache (refetch)", variable=self.refresh_var
        ).grid(row=len(actions), column=0, sticky="w", pady=(12, 0))

        ttk.Button(controls, text="Export result JSON", width=22, command=self.export_json).grid(
            row=len(actions) + 1, column=0, sticky="ew", pady=(12, 2)
        )

        self._outcome: StructureSearchOutcome | None = None
        self._local_pdb: str | None = None

    # -- actions ------------------------------------------------------------

    def open_fasta(self) -> None:
        path = filedialog.askopenfilename(
            filetypes=[("FASTA", "*.fasta *.fa *.faa *.txt"), ("All files", "*.*")]
        )
        if path:
            with open(path, encoding="utf-8", errors="replace") as handle:
                self.input.set(handle.read())

    def set_query(self, text: str) -> None:
        self.input.set(text)

    def run_retrieve(self, providers: list[str] | None = None, stop_on_exact: bool = True) -> None:
        raw = self.input.get().strip()
        if not raw:
            messagebox.showinfo("BioSeqInsight", "Enter a sequence, accession or PDB identifier.")
            return
        records = parse_fasta(raw)
        query = records[0].sequence if records and records[0].sequence else raw
        refresh = self.refresh_var.get()

        self.app.status.busy("Querying structural resources...")
        self.output.set("Working. External services can take up to a minute.\n")

        def work():
            return self.app.manager.retrieve(
                query, providers=providers, stop_on_exact=stop_on_exact, refresh=refresh
            )

        def done(outcome: StructureSearchOutcome) -> None:
            self._outcome = outcome
            self.output.set(format_structure_outcome(outcome))
            self.badge.show(outcome.best_mapping.value, outcome.best_mapping.label)
            best = outcome.best
            if best and best.download_path:
                self._local_pdb = best.download_path
            self.app.status.idle(
                f"Finished in {outcome.elapsed_s:.1f}s - best match {outcome.best_mapping.value}."
            )

        def failed(exc: Exception) -> None:
            self.output.set(f"Retrieval failed.\n\n{type(exc).__name__}: {exc}")
            self.app.status.idle("Retrieval failed.")

        self.runner.submit(work, done, failed)

    def run_retrieve_all(self) -> None:
        self.run_retrieve(providers=None, stop_on_exact=False)

    def load_local(self) -> None:
        path = filedialog.askopenfilename(filetypes=[("PDB", "*.pdb *.ent"), ("All files", "*.*")])
        if not path:
            return
        self._local_pdb = path
        text = Path(path).read_text(encoding="utf-8", errors="replace")
        from ..structures.pdbio import parse_pdb

        parsed = parse_pdb(text)
        lines = [
            "Local coordinate file",
            DIVIDER,
            f"Path          {path}",
            f"Atoms         {parsed.atom_count}",
            f"Chains        {len(parsed.chains)}",
        ]
        for chain in parsed.chains:
            lines.append(f"  chain {chain.chain_id}: {chain.residue_count} residues ({chain.source})")
        if parsed.accessions:
            lines.append(f"Cross-refs    {', '.join(parsed.accessions)}")
        if parsed.title:
            lines.append(f"Title         {parsed.title}")
        lines += [
            "",
            "A locally loaded file is not identity-checked against any query.",
        ]
        self.output.set("\n".join(lines))
        self.badge.show("M1", "Local file, identity unverified")

    def view_browser(self) -> None:
        path = self._current_pdb()
        if not path:
            return
        text = Path(path).read_text(encoding="utf-8", errors="replace")
        best = self._outcome.best if self._outcome else None
        html = write_viewer_html(
            text,
            title=Path(path).name,
            source=best.source if best else "local file",
            confidence=(
                f"{best.confidence_kind} {best.confidence:.1f}"
                if best and best.confidence is not None
                else ""
            ),
            mapping=(
                f"{best.mapping_level.value} {best.mapping_level.label}" if best else ""
            ),
        )
        if open_in_browser(html):
            self.app.status.idle("Opened the 3Dmol.js viewer in your browser.")
        else:
            messagebox.showinfo("BioSeqInsight", f"Viewer page written to {html}")

    def view_external(self) -> None:
        path = self._current_pdb()
        if path:
            self.app.status.idle(
                open_in_external_viewer(path, preferred=self.app.settings.viewer_command)
            )

    def _current_pdb(self) -> str | None:
        if self._local_pdb and Path(self._local_pdb).is_file():
            return self._local_pdb
        messagebox.showinfo(
            "BioSeqInsight", "Retrieve a structure or load a local PDB file first."
        )
        return None

    def export_json(self) -> None:
        if self._outcome is None:
            messagebox.showinfo("BioSeqInsight", "Retrieve a structure first.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".json")
        if path:
            Path(path).write_text(to_json(self._outcome), encoding="utf-8")
            self.app.status.idle(f"Exported to {path}")


class BatchView(ttk.Frame):
    """The Batch Analysis tab."""

    def __init__(self, master, app, **kwargs):
        super().__init__(master, padding=12, **kwargs)
        self.app = app
        self.runner = BackgroundRunner(self)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(2, weight=1)

        top = ttk.Frame(self)
        top.grid(row=0, column=0, sticky="ew")
        top.columnconfigure(1, weight=1)
        ttk.Label(top, text="FASTA file").grid(row=0, column=0, sticky="w")
        self.path_var = tk.StringVar()
        ttk.Entry(top, textvariable=self.path_var).grid(row=0, column=1, sticky="ew", padx=8)
        ttk.Button(top, text="Browse...", command=self.browse).grid(row=0, column=2)

        options = ttk.LabelFrame(self, text="Options", padding=10)
        options.grid(row=1, column=0, sticky="ew", pady=10)
        self.structures_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            options,
            text="Retrieve and validate structures (requires network)",
            variable=self.structures_var,
        ).grid(row=0, column=0, columnspan=4, sticky="w")
        ttk.Label(options, text="Workers").grid(row=1, column=0, sticky="w", pady=(8, 0))
        self.workers = tk.IntVar(value=4)
        ttk.Spinbox(options, from_=1, to=16, textvariable=self.workers, width=6).grid(
            row=1, column=1, sticky="w", pady=(8, 0)
        )
        ttk.Label(options, text="Min ORF (aa)").grid(row=1, column=2, sticky="e", pady=(8, 0), padx=(16, 4))
        self.min_orf = tk.IntVar(value=30)
        ttk.Spinbox(options, from_=1, to=2000, textvariable=self.min_orf, width=6).grid(
            row=1, column=3, sticky="w", pady=(8, 0)
        )
        ttk.Button(options, text="Run batch", command=self.run).grid(
            row=2, column=0, sticky="w", pady=(12, 0)
        )
        for index, (label, kind) in enumerate(
            [("Export CSV", "csv"), ("Export TSV", "tsv"), ("Export JSON", "json"), ("Export HTML", "html")]
        ):
            ttk.Button(options, text=label, command=lambda k=kind: self.export(k)).grid(
                row=2, column=index + 1 if index < 3 else 1, sticky="w", pady=(12, 0), padx=4
            )

        self.output = LabelledText(self, "Batch log", height=18, readonly=True)
        self.output.grid(row=2, column=0, sticky="nsew")
        bind_copy(self.output.text)

        self._rows = []
        self._summary = None

    def browse(self) -> None:
        path = filedialog.askopenfilename(
            filetypes=[("FASTA", "*.fasta *.fa *.faa *.fna *.txt"), ("All files", "*.*")]
        )
        if path:
            self.path_var.set(path)

    def run(self) -> None:
        path = self.path_var.get().strip()
        if not path or not Path(path).is_file():
            messagebox.showinfo("BioSeqInsight", "Choose a FASTA file first.")
            return
        options = BatchOptions(
            include_structures=self.structures_var.get(),
            workers=self.workers.get(),
            min_orf_aa=self.min_orf.get(),
        )
        self.output.set(f"Reading {path}\n")
        self.app.status.busy("Running batch analysis...")

        def work():
            return self.app.batch_runner.run_file(path, options)

        def done(payload) -> None:
            rows, summary = payload
            self._rows, self._summary = rows, summary
            self.output.append(self._format_summary(rows, summary))
            self.app.status.idle(f"Batch finished: {summary.total} record(s).")

        def failed(exc: Exception) -> None:
            self.output.append(f"\nBatch failed: {type(exc).__name__}: {exc}\n")
            self.app.status.idle("Batch failed.")

        self.runner.submit(work, done, failed)

    @staticmethod
    def _format_summary(rows, summary) -> str:
        lines = [
            "",
            "Batch summary",
            DIVIDER,
            f"Records          {summary.total}",
            f"Succeeded        {summary.succeeded}",
            f"Failed           {summary.failed}",
            f"Elapsed          {summary.elapsed_s:.2f} s",
        ]
        if summary.elapsed_s > 0:
            lines.append(f"Throughput       {summary.total / summary.elapsed_s:.1f} records/s")
        if summary.structures_attempted:
            lines += [
                f"Structures tried {summary.structures_attempted}",
                f"Structures found {summary.structures_retrieved}",
                f"Cache hits       {summary.cache_hits}",
                f"Recovered calls  {summary.recovered_requests}",
            ]
        if summary.mapping_counts:
            lines.append("Mapping levels   " + ", ".join(
                f"{level}={count}" for level, count in sorted(summary.mapping_counts.items())
            ))
        failures = [row for row in rows if row.error]
        if failures:
            lines += ["", "Records with problems:"]
            lines += [f"  {row.identifier}: {row.error}" for row in failures[:20]]
        return "\n".join(lines) + "\n"

    def export(self, kind: str) -> None:
        if not self._rows:
            messagebox.showinfo("BioSeqInsight", "Run a batch before exporting.")
            return
        path = filedialog.asksaveasfilename(defaultextension=f".{kind}")
        if not path:
            return
        if kind == "csv":
            write_csv(self._rows, path)
        elif kind == "tsv":
            write_tsv(self._rows, path)
        elif kind == "html":
            write_html_report(self._rows, path)
        else:
            write_json({"rows": [row.to_dict() for row in self._rows]}, path)
        self.app.status.idle(f"Exported to {path}")


__all__ = ["BatchView", "StructureView", "format_structure_outcome"]
