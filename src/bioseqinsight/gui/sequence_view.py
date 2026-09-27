"""Sequence and protein tabs.

Both views are thin: they collect text, call the core API, and hand the
returned dataclass to :mod:`bioseqinsight.gui.formatting`. No calculation
happens here, and the formatting functions live in a Tkinter-free module so
that the displayed wording is tested headless.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from ..core.alphabet import SequenceError, detect_alphabet, parse_fasta
from ..core.protein import ProteinAnalyzer, hydropathy_profile
from ..core.sequence import SequenceAnalyzer
from ..core.translation import six_frame_translation
from ..io.export import (
    protein_result_to_record,
    sequence_result_to_record,
    to_json,
    write_records_csv,
)
from ..models.results import ProteinResult, SequenceResult
from .formatting import (
    DIVIDER,
    _wrap,
    format_hydropathy,
    format_protein_result,
    format_sequence_result,
    format_six_frames,
)
from .widgets import LabelledText, bind_copy

# --------------------------------------------------------------------------
# views
# --------------------------------------------------------------------------


class SequenceView(ttk.Frame):
    """The Sequence Analysis tab."""

    def __init__(self, master, app, **kwargs):
        super().__init__(master, padding=12, **kwargs)
        self.app = app
        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)

        self.input = LabelledText(self, "Input DNA / RNA sequence or FASTA", height=7)
        self.input.grid(row=0, column=0, columnspan=2, sticky="nsew", pady=(0, 10))

        controls = ttk.Frame(self)
        controls.grid(row=1, column=0, sticky="nw", padx=(0, 12))
        self.output = LabelledText(self, "Result", height=18, readonly=True)
        self.output.grid(row=1, column=1, sticky="nsew")
        bind_copy(self.output.text)

        buttons = [
            ("Open FASTA...", self.open_fasta),
            ("Full analysis", self.run_full),
            ("GC content", self.run_gc),
            ("Reverse complement", self.run_revcomp),
            ("Transcribe to mRNA", self.run_transcribe),
            ("Six-frame translation", self.run_six_frames),
            ("Find ORFs", self.run_orfs),
            ("Melting temperature", self.run_tm),
        ]
        for index, (label, command) in enumerate(buttons):
            ttk.Button(controls, text=label, width=22, command=command).grid(
                row=index, column=0, sticky="ew", pady=2
            )

        ttk.Label(controls, text="Motif (IUPAC)").grid(row=len(buttons), column=0, sticky="w", pady=(12, 2))
        self.motif_var = tk.StringVar()
        ttk.Entry(controls, textvariable=self.motif_var, width=24).grid(
            row=len(buttons) + 1, column=0, sticky="ew"
        )
        ttk.Button(controls, text="Find motif", width=22, command=self.run_motif).grid(
            row=len(buttons) + 2, column=0, sticky="ew", pady=2
        )

        ttk.Label(controls, text="Minimum ORF length (aa)").grid(
            row=len(buttons) + 3, column=0, sticky="w", pady=(12, 2)
        )
        self.min_orf = tk.IntVar(value=30)
        ttk.Spinbox(controls, from_=1, to=2000, textvariable=self.min_orf, width=8).grid(
            row=len(buttons) + 4, column=0, sticky="w"
        )

        export = ttk.Frame(controls)
        export.grid(row=len(buttons) + 5, column=0, sticky="ew", pady=(16, 0))
        ttk.Button(export, text="Export CSV", width=10, command=lambda: self.export("csv")).grid(
            row=0, column=0, padx=(0, 4)
        )
        ttk.Button(export, text="Export JSON", width=11, command=lambda: self.export("json")).grid(
            row=0, column=1
        )

        self._last: SequenceResult | None = None

    # -- helpers ------------------------------------------------------------

    def _analyzer(self) -> SequenceAnalyzer:
        text = self.input.get().strip()
        if not text:
            raise SequenceError("Paste a sequence or open a FASTA file first.")
        records = parse_fasta(text)
        identifier = records[0].identifier if records else "sequence"
        sequence = records[0].sequence if records else text
        if detect_alphabet(sequence) == "protein":
            raise SequenceError(
                "That looks like a protein sequence. Use the Protein Analysis tab."
            )
        return SequenceAnalyzer(
            sequence, identifier=identifier, table=self.app.settings.genetic_code_table
        )

    def _guard(self, function):
        try:
            function()
        except SequenceError as exc:
            self.output.set(str(exc))
            self.app.status.idle("Input rejected.")
        except Exception as exc:  # pragma: no cover - defensive
            messagebox.showerror("BioSeqInsight", f"{type(exc).__name__}: {exc}")
            self.app.status.idle("Error.")

    # -- actions ------------------------------------------------------------

    def open_fasta(self) -> None:
        path = filedialog.askopenfilename(
            title="Open FASTA",
            filetypes=[("FASTA", "*.fasta *.fa *.fna *.txt"), ("All files", "*.*")],
        )
        if not path:
            return
        with open(path, encoding="utf-8", errors="replace") as handle:
            self.input.set(handle.read())
        self.app.status.idle(f"Loaded {path}")

    def run_full(self) -> None:
        def work():
            analyzer = self._analyzer()
            result = analyzer.analyze(
                motif=self.motif_var.get().strip() or None, min_orf_aa=self.min_orf.get()
            )
            self._last = result
            self.output.set(format_sequence_result(result))
            self.app.status.idle(f"Analysed {result.length} nt.")

        self._guard(work)

    def run_gc(self) -> None:
        def work():
            analyzer = self._analyzer()
            composition = analyzer.get_composition()
            self.output.set(
                "\n".join(
                    [
                        "GC content",
                        DIVIDER,
                        f"GC content   {analyzer.get_gc_content():.2f}%",
                        f"Length       {len(analyzer)} nt",
                        "Counts       "
                        + ", ".join(f"{b}={c}" for b, c in composition.items()),
                        "",
                        "Ambiguity codes are excluded from the GC calculation.",
                    ]
                )
            )
            self.app.status.idle("GC content computed.")

        self._guard(work)

    def run_revcomp(self) -> None:
        self._guard(
            lambda: self.output.set(
                "Reverse complement\n" + DIVIDER + "\n" + _wrap(self._analyzer().get_reverse_complement(), 60)
            )
        )

    def run_transcribe(self) -> None:
        self._guard(
            lambda: self.output.set(
                "mRNA (coding strand)\n" + DIVIDER + "\n" + _wrap(self._analyzer().transcribe(), 60)
            )
        )

    def run_six_frames(self) -> None:
        self._guard(
            lambda: self.output.set(
                format_six_frames(
                    six_frame_translation(
                        self._analyzer().sequence, table=self.app.settings.genetic_code_table
                    )
                )
            )
        )

    def run_orfs(self) -> None:
        def work():
            analyzer = self._analyzer()
            result = analyzer.analyze(min_orf_aa=self.min_orf.get(), include_sequences=False)
            self._last = result
            self.output.set(format_sequence_result(result))

        self._guard(work)

    def run_tm(self) -> None:
        def work():
            tm = self._analyzer().calculate_tm(
                primer_nM=self.app.settings.tm_primer_nM, na_mM=self.app.settings.tm_na_mM
            )
            lines = [
                "Melting temperature",
                DIVIDER,
                f"Tm       {tm.tm_c:.2f} C",
                f"Method   {tm.method}",
                f"Length   {tm.length} nt",
            ]
            if tm.method == "nearest_neighbor":
                lines.append(f"Wallace value would be {tm.wallace_c:.0f} C and is not applicable.")
            if tm.note:
                lines += ["", tm.note]
            self.output.set("\n".join(lines))

        self._guard(work)

    def run_motif(self) -> None:
        def work():
            motif = self.motif_var.get().strip()
            if not motif:
                raise SequenceError("Enter a motif, for example GAATTC or WGATAR.")
            hits = self._analyzer().find_motifs(motif)
            lines = [f"Motif search: {motif}", DIVIDER, f"{len(hits)} hit(s)", ""]
            lines += [f"  {h.strand} strand  {h.start}-{h.end}" for h in hits[:200]]
            lines += ["", "Exact IUPAC matching, not statistical motif discovery."]
            self.output.set("\n".join(lines))

        self._guard(work)

    def export(self, kind: str) -> None:
        if self._last is None:
            messagebox.showinfo("BioSeqInsight", "Run an analysis before exporting.")
            return
        extension = ".csv" if kind == "csv" else ".json"
        path = filedialog.asksaveasfilename(defaultextension=extension)
        if not path:
            return
        if kind == "csv":
            write_records_csv([sequence_result_to_record(self._last)], path)
        else:
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(to_json(self._last))
        self.app.status.idle(f"Exported to {path}")


class ProteinView(ttk.Frame):
    """The Protein Analysis tab."""

    def __init__(self, master, app, **kwargs):
        super().__init__(master, padding=12, **kwargs)
        self.app = app
        self.columnconfigure(1, weight=1)
        self.rowconfigure(1, weight=1)

        self.input = LabelledText(self, "Input protein sequence or FASTA", height=7)
        self.input.grid(row=0, column=0, columnspan=2, sticky="nsew", pady=(0, 10))

        controls = ttk.Frame(self)
        controls.grid(row=1, column=0, sticky="nw", padx=(0, 12))
        self.output = LabelledText(self, "Result", height=18, readonly=True)
        self.output.grid(row=1, column=1, sticky="nsew")
        bind_copy(self.output.text)

        buttons = [
            ("Open FASTA...", self.open_fasta),
            ("Full analysis", self.run_full),
            ("Composition", self.run_composition),
            ("Hydropathy profile", self.run_hydropathy),
            ("Propensity sketch", self.run_sketch),
            ("Send to Structure tab", self.send_to_structure),
        ]
        for index, (label, command) in enumerate(buttons):
            ttk.Button(controls, text=label, width=22, command=command).grid(
                row=index, column=0, sticky="ew", pady=2
            )

        ttk.Label(controls, text="Hydropathy window").grid(
            row=len(buttons), column=0, sticky="w", pady=(12, 2)
        )
        self.window = tk.IntVar(value=9)
        ttk.Spinbox(controls, from_=3, to=51, increment=2, textvariable=self.window, width=8).grid(
            row=len(buttons) + 1, column=0, sticky="w"
        )

        export = ttk.Frame(controls)
        export.grid(row=len(buttons) + 2, column=0, sticky="ew", pady=(16, 0))
        ttk.Button(export, text="Export CSV", width=10, command=lambda: self.export("csv")).grid(
            row=0, column=0, padx=(0, 4)
        )
        ttk.Button(export, text="Export JSON", width=11, command=lambda: self.export("json")).grid(
            row=0, column=1
        )

        self._last: ProteinResult | None = None

    def _analyzer(self) -> ProteinAnalyzer:
        text = self.input.get().strip()
        if not text:
            raise SequenceError("Paste a protein sequence or open a FASTA file first.")
        records = parse_fasta(text)
        identifier = records[0].identifier if records else "protein"
        sequence = records[0].sequence if records else text
        if detect_alphabet(sequence) == "dna":
            raise SequenceError(
                "That looks like a nucleotide sequence. Translate it on the "
                "Sequence Analysis tab first."
            )
        return ProteinAnalyzer(sequence, identifier=identifier)

    def _guard(self, function):
        try:
            function()
        except SequenceError as exc:
            self.output.set(str(exc))
        except Exception as exc:  # pragma: no cover - defensive
            messagebox.showerror("BioSeqInsight", f"{type(exc).__name__}: {exc}")

    def open_fasta(self) -> None:
        path = filedialog.askopenfilename(
            title="Open FASTA", filetypes=[("FASTA", "*.fasta *.fa *.faa *.txt"), ("All files", "*.*")]
        )
        if not path:
            return
        with open(path, encoding="utf-8", errors="replace") as handle:
            self.input.set(handle.read())

    def run_full(self) -> None:
        def work():
            result = self._analyzer().analyze(include_sketch=True)
            self._last = result
            self.output.set(format_protein_result(result))
            self.app.status.idle(f"Analysed {result.length} aa.")

        self._guard(work)

    def run_composition(self) -> None:
        def work():
            result = self._analyzer().analyze(include_sketch=False)
            self._last = result
            self.output.set(format_protein_result(result))

        self._guard(work)

    def run_hydropathy(self) -> None:
        self._guard(
            lambda: self.output.set(
                format_hydropathy(
                    hydropathy_profile(self._analyzer().sequence, window=self.window.get()),
                    self.window.get(),
                )
            )
        )

    def run_sketch(self) -> None:
        def work():
            sketch = self._analyzer().get_secondary_structure_propensity()
            self.output.set(
                "\n".join(
                    [
                        "Secondary-structure propensity sketch",
                        DIVIDER,
                        sketch.method,
                        f"helix {sketch.helix_percent:.1f}%   sheet {sketch.sheet_percent:.1f}%   "
                        f"coil {sketch.coil_percent:.1f}%",
                        "",
                        _wrap(sketch.visual, 60),
                        "",
                        "This is an illustrative sketch. For a real assignment use DSSP on a "
                        "structure, or a dedicated predictor.",
                    ]
                )
            )

        self._guard(work)

    def send_to_structure(self) -> None:
        self._guard(lambda: self.app.send_to_structure(self._analyzer().sequence))

    def export(self, kind: str) -> None:
        if self._last is None:
            messagebox.showinfo("BioSeqInsight", "Run an analysis before exporting.")
            return
        extension = ".csv" if kind == "csv" else ".json"
        path = filedialog.asksaveasfilename(defaultextension=extension)
        if not path:
            return
        if kind == "csv":
            write_records_csv([protein_result_to_record(self._last)], path)
        else:
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(to_json(self._last))
        self.app.status.idle(f"Exported to {path}")


__all__ = [
    "ProteinView",
    "SequenceView",
    "format_hydropathy",
    "format_protein_result",
    "format_sequence_result",
    "format_six_frames",
]
