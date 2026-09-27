"""The BioSeqInsight main window.

Responsibilities are deliberately narrow: build the window, own the shared
service objects (settings, structure manager, batch runner, current project),
and route menu commands to the views. Every biological computation and every
network request happens in the layers underneath.
"""

from __future__ import annotations

import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

from .. import __version__
from ..config.settings import Settings, get_settings
from ..services.logging_setup import configure_logging
from ..structures.manager import StructureManager
from ..workflows.batch import BatchRunner
from ..workflows.project import Project, ProjectError
from .sequence_view import ProteinView, SequenceView
from .structure_view import BatchView, StructureView
from .widgets import StatusBar

APP_TITLE = "BioSeqInsight"


class BioSeqInsightApp(ttk.Frame):
    """Root application frame."""

    def __init__(self, master: tk.Tk, settings: Settings | None = None):
        super().__init__(master)
        self.master = master
        self.settings = settings or get_settings()
        self.settings.ensure_directories()
        configure_logging(
            log_file=self.settings.log_file,
            level=self.settings.log_level,
            json_format=self.settings.log_json,
            console=False,
        )

        self.manager = StructureManager(self.settings)
        self.batch_runner = BatchRunner(self.settings, manager=self.manager)
        self.project: Project | None = None

        master.title(f"{APP_TITLE} {__version__}")
        master.geometry("1180x760")
        master.minsize(900, 600)
        self._style()

        self.pack(fill="both", expand=True)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.notebook = ttk.Notebook(self)
        self.notebook.grid(row=0, column=0, sticky="nsew")
        self.sequence_view = SequenceView(self.notebook, self)
        self.protein_view = ProteinView(self.notebook, self)
        self.structure_view = StructureView(self.notebook, self)
        self.batch_view = BatchView(self.notebook, self)
        self.notebook.add(self.sequence_view, text="Sequence Analysis")
        self.notebook.add(self.protein_view, text="Protein Analysis")
        self.notebook.add(self.structure_view, text="Structure Retrieval")
        self.notebook.add(self.batch_view, text="Batch Analysis")

        self.status = StatusBar(self, padding=(10, 6))
        self.status.grid(row=1, column=0, sticky="ew")
        self._build_menu()
        self.status.idle(
            f"Ready. Settings from {self.settings.source}. "
            + ("Offline mode." if self.settings.offline else "Network enabled.")
        )

    # -- chrome -------------------------------------------------------------

    def _style(self) -> None:
        style = ttk.Style()
        if "clam" in style.theme_names():
            style.theme_use("clam")
        style.configure("TNotebook.Tab", padding=(16, 8))
        style.configure("TButton", padding=(8, 4))

    def _build_menu(self) -> None:
        menubar = tk.Menu(self.master)

        project_menu = tk.Menu(menubar, tearoff=0)
        project_menu.add_command(label="New project...", command=self.new_project)
        project_menu.add_command(label="Open project...", command=self.open_project)
        project_menu.add_command(label="Export project as zip...", command=self.export_project)
        project_menu.add_separator()
        project_menu.add_command(label="Project info", command=self.show_project)
        project_menu.add_separator()
        project_menu.add_command(label="Quit", command=self.master.destroy)
        menubar.add_cascade(label="Project", menu=project_menu)

        tools = tk.Menu(menubar, tearoff=0)
        tools.add_command(label="Cache statistics", command=self.show_cache)
        tools.add_command(label="Clear cache", command=self.clear_cache)
        tools.add_separator()
        tools.add_command(label="Settings...", command=self.show_settings)
        menubar.add_cascade(label="Tools", menu=tools)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="Mapping levels (M0-M4)", command=self.show_mapping_help)
        help_menu.add_command(label="About BioSeqInsight", command=self.show_about)
        menubar.add_cascade(label="Help", menu=help_menu)

        self.master.config(menu=menubar)

    # -- cross-tab wiring ---------------------------------------------------

    def send_to_structure(self, sequence: str) -> None:
        self.structure_view.set_query(sequence)
        self.notebook.select(self.structure_view)
        self.status.idle("Sequence copied to the Structure Retrieval tab.")

    # -- project commands ---------------------------------------------------

    def new_project(self) -> None:
        path = filedialog.askdirectory(title="Choose an empty folder for the new project")
        if not path:
            return
        try:
            self.project = Project.create(path, name=Path(path).name, settings=self.settings, overwrite=True)
        except ProjectError as exc:
            messagebox.showerror(APP_TITLE, str(exc))
            return
        self.status.idle(f"Project created at {path}")

    def open_project(self) -> None:
        path = filedialog.askdirectory(title="Open a BioSeqInsight project folder")
        if not path:
            return
        try:
            self.project = Project.open(path)
        except ProjectError as exc:
            messagebox.showerror(APP_TITLE, str(exc))
            return
        self.status.idle(f"Opened project {self.project.manifest.name}")

    def export_project(self) -> None:
        if self.project is None:
            messagebox.showinfo(APP_TITLE, "Create or open a project first.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".zip")
        if path:
            messagebox.showinfo(APP_TITLE, f"Project exported to\n{self.project.export_zip(path)}")

    def show_project(self) -> None:
        if self.project is None:
            messagebox.showinfo(APP_TITLE, "No project is open.")
            return
        summary = self.project.summary()
        messagebox.showinfo(
            APP_TITLE, "\n".join(f"{key}: {value}" for key, value in summary.items())
        )

    # -- tools --------------------------------------------------------------

    def show_cache(self) -> None:
        stats = self.manager.cache_stats()
        messagebox.showinfo(
            APP_TITLE, "\n".join(f"{key}: {value}" for key, value in stats.items())
        )

    def clear_cache(self) -> None:
        if messagebox.askyesno(APP_TITLE, "Delete every cached response?"):
            removed = self.manager.cache.clear()
            self.status.idle(f"Cleared {removed} cached payload(s).")

    def show_settings(self) -> None:
        window = tk.Toplevel(self.master)
        window.title("Settings")
        window.geometry("560x420")
        text = tk.Text(window, wrap="none")
        text.pack(fill="both", expand=True)
        import json

        text.insert(
            "1.0",
            json.dumps({"source": self.settings.source, **self.settings.to_dict()}, indent=2, sort_keys=True),
        )
        text.configure(state="disabled")
        ttk.Label(
            window,
            text=(
                "Settings are read from ~/.bioseqinsight/settings.json, environment "
                "variables and command-line flags. Edit the file to change them."
            ),
            wraplength=540,
            padding=8,
        ).pack(fill="x")

    # -- help ---------------------------------------------------------------

    def show_mapping_help(self) -> None:
        messagebox.showinfo(
            "Sequence-to-structure mapping levels",
            "M4  Exact: the queried accession is cross-referenced by the structure and\n"
            "    sequence identity and coverage meet the exact thresholds.\n\n"
            "M3  High confidence: identity and coverage are high but the accession was\n"
            "    not confirmed.\n\n"
            "M2  Related record only: a homologue, orthologue or fragment. This is NOT\n"
            "    the structure of your protein.\n\n"
            "M1  Retrieved but unverified: no comparison was possible.\n\n"
            "M0  No structural record was returned.\n\n"
            "Only M4 justifies stating that a retrieved structure is the requested "
            "molecule.",
        )

    def show_about(self) -> None:
        messagebox.showinfo(
            f"About {APP_TITLE}",
            f"BioSeqInsight {__version__}\n\n"
            "An integrated desktop workflow for DNA sequence analysis and validated "
            "protein structure-resource retrieval.\n\n"
            "Structures come from RCSB PDB, the AlphaFold Protein Structure Database "
            "and the ESM Atlas fold service. BioSeqInsight does not implement or "
            "redistribute AlphaFold or ESMFold; please cite those resources when you "
            "use structures obtained through this software.\n\n"
            "MIT licence.",
        )


def launch(settings: Settings | None = None) -> int:
    """Start the GUI. Returns a process exit code."""
    try:
        root = tk.Tk()
    except tk.TclError as exc:  # pragma: no cover - headless environment
        print(
            f"Cannot start the graphical interface: {exc}\n"
            "Use the command-line interface instead: python -m bioseqinsight --help",
        )
        return 1
    BioSeqInsightApp(root, settings=settings)
    root.mainloop()
    return 0


def main() -> int:  # pragma: no cover - entry point
    return launch()


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
