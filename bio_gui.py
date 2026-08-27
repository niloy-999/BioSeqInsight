#!/usr/bin/env python3
"""Tkinter front-end for BioSeqInsight."""

from __future__ import annotations

import os
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, scrolledtext, ttk

from Bio import SeqIO
from PIL import Image, ImageTk

import sequence_operations as seq_ops
import structure_predictions as pred_ops


APP_TITLE = "BioSeqInsight"
WINDOW_SIZE = "1000x850"


def _read_first_fasta_record(kind: str) -> str:
    title = f"Select {kind} FASTA file"
    path = filedialog.askopenfilename(
        title=title,
        filetypes=(("FASTA files", "*.fasta *.fa *.fna *.faa"), ("All files", "*.*")),
    )
    if not path:
        return ""
    try:
        record = next(SeqIO.parse(path, "fasta"))
    except Exception as exc:
        messagebox.showerror("FASTA error", str(exc))
        return ""
    return str(record.seq)


def _background_path() -> str | None:
    names = ("background.jpg", "background.png", "background.jpeg")
    folders = (Path.cwd(), Path(__file__).resolve().parent)
    for folder in folders:
        for name in names:
            candidate = folder / name
            if candidate.is_file():
                return str(candidate)
    return None


def create_main_window() -> None:
    root = tk.Tk()
    root.title(APP_TITLE)
    root.geometry(WINDOW_SIZE)
    root.minsize(900, 640)

    notebook = ttk.Notebook(root)
    notebook.pack(expand=True, fill="both")

    # tk.Frame, not ttk.Frame: themed frames paint an opaque panel over the photo.
    seq_tab = tk.Frame(notebook)
    pred_tab = tk.Frame(notebook)
    notebook.add(seq_tab, text="Sequence Analysis")
    notebook.add(pred_tab, text="3D Structure Prediction")

    add_background_image(seq_tab)
    add_background_image(pred_tab)
    add_sequence_analysis_tab(seq_tab)
    add_prediction_tab(pred_tab)

    # Keep the photo behind widgets after they are packed.
    seq_tab.after(100, lambda: _lower_background(seq_tab))
    pred_tab.after(100, lambda: _lower_background(pred_tab))

    root.mainloop()


def _lower_background(tab: tk.Frame) -> None:
    label = getattr(tab, "_bg_label", None)
    if label is not None:
        label.lower()


def add_background_image(tab: tk.Frame) -> None:
    path = _background_path()
    if path is None:
        print(
            "Background image not found. Put background.jpg in the same folder "
            "as bio_gui.py or in the folder you launch Python from."
        )
        return

    tab._bg_path = path
    tab._bg_label = tk.Label(tab, borderwidth=0, highlightthickness=0)
    tab._bg_label.place(x=0, y=0, relwidth=1, relheight=1)
    tab._bg_job = None

    def refresh(_event=None) -> None:
        width = max(tab.winfo_width(), 1)
        height = max(tab.winfo_height(), 1)
        if width < 20 or height < 20:
            return
        if tab._bg_job is not None:
            tab.after_cancel(tab._bg_job)

        def paint() -> None:
            image = Image.open(tab._bg_path).resize((width, height), Image.LANCZOS)
            photo = ImageTk.PhotoImage(image)
            tab._bg_photo = photo
            tab._bg_label.configure(image=photo)
            tab._bg_label.lower()

        tab._bg_job = tab.after(40, paint)

    tab.bind("<Configure>", refresh)
    tab.after(80, refresh)


def _make_io_boxes(parent: tk.Frame, input_caption: str):
    """Pack controls on the tab itself so the wallpaper stays visible."""
    tk.Label(parent, text=input_caption).pack(anchor="w", padx=5, pady=(6, 0))
    input_text = scrolledtext.ScrolledText(parent, height=4, wrap=tk.WORD)
    input_text.pack(fill="x", padx=5, pady=5)

    buttons = tk.Frame(parent)
    buttons.pack(anchor="w", padx=5, pady=4)

    tk.Label(parent, text="Output:").pack(anchor="w", padx=5)
    output_text = scrolledtext.ScrolledText(parent, height=6, wrap=tk.WORD, state="disabled")
    output_text.pack(fill="x", padx=5, pady=5)

    def show(title: str, result: str) -> None:
        output_text.config(state="normal")
        output_text.delete("1.0", tk.END)
        output_text.insert(tk.END, f"{title}\n{'-' * len(title)}\n{result}\n")
        output_text.config(state="disabled")

    def raw_input() -> str:
        return input_text.get("1.0", tk.END)

    return input_text, buttons, show, raw_input


def add_sequence_analysis_tab(tab: tk.Frame) -> None:
    input_text, buttons, show, raw_input = _make_io_boxes(tab, "Input DNA Sequence:")

    def upload() -> None:
        seq = _read_first_fasta_record("DNA")
        if seq:
            input_text.delete("1.0", tk.END)
            input_text.insert(tk.END, seq)

    rows = [
        ("Upload DNA FASTA File", upload),
        ("Calculate GC Content", lambda: show("GC Content", seq_ops.format_gc_report(raw_input()))),
        ("Reverse Sequence", lambda: show("Reversed Sequence", seq_ops.reverse_sequence(raw_input()) or "No DNA found.")),
        (
            "Reverse Complement",
            lambda: show("Reverse Complement", seq_ops.reverse_complement(raw_input()) or "No DNA found."),
        ),
    ]
    for text, cmd in rows:
        tk.Button(buttons, text=text, command=cmd).pack(anchor="w", padx=5, pady=2)

    tk.Label(buttons, text="Motif:").pack(anchor="w", padx=5, pady=(6, 0))
    motif_entry = tk.Entry(buttons)
    motif_entry.pack(anchor="w", padx=5, pady=2)
    tk.Button(
        buttons,
        text="Find Motif",
        command=lambda: show("Motif Positions", seq_ops.format_motif_report(raw_input(), motif_entry.get())),
    ).pack(anchor="w", padx=5, pady=2)

    more = [
        ("Nucleotide Frequency", lambda: show("Nucleotide Frequency", seq_ops.format_frequency_report(raw_input()))),
        ("Melting Temperature (Tm)", lambda: show("Melting Temperature", seq_ops.format_tm_report(raw_input()))),
        ("Transcript to mRNA", lambda: show("mRNA Sequence", seq_ops.dna_to_mrna(raw_input()) or "No DNA found.")),
        ("Translate to Protein", lambda: show("Protein Sequence", seq_ops.format_translate_report(raw_input()))),
        ("Translate with ORF", lambda: show("Open Reading Frames", seq_ops.format_orf_report(raw_input()))),
    ]
    for text, cmd in more:
        tk.Button(buttons, text=text, command=cmd).pack(anchor="w", padx=5, pady=2)


def add_prediction_tab(tab: tk.Frame) -> None:
    input_text, buttons, show, raw_input = _make_io_boxes(tab, "Input Protein Sequence for Analysis:")
    status = tk.Label(tab, text="Ready.")
    status.pack(anchor="w", padx=5, pady=(0, 8))

    def upload() -> None:
        seq = _read_first_fasta_record("protein")
        if seq:
            input_text.delete("1.0", tk.END)
            input_text.insert(tk.END, seq)

    def set_busy(busy: bool, text: str) -> None:
        status.config(text=text)
        state = "disabled" if busy else "normal"
        for child in buttons.winfo_children():
            try:
                child.config(state=state)
            except tk.TclError:
                pass

    def run_job(label: str, fn) -> None:
        set_busy(True, label)

        def work():
            result = fn()
            buttons.after(0, lambda: finish(result))

        def finish(result: str) -> None:
            set_busy(False, "Ready.")
            show(label, result)

        threading.Thread(target=work, daemon=True).start()

    def load_pdb() -> None:
        path = filedialog.askopenfilename(
            title="Select PDB file",
            filetypes=(("PDB files", "*.pdb *.ent"), ("All files", "*.*")),
        )
        if path:
            show("Load PDB", pred_ops.load_pdb_file(path))

    tk.Button(buttons, text="Upload Protein FASTA File", command=upload).pack(anchor="w", padx=5, pady=2)
    tk.Button(
        buttons,
        text="Predict 3D Structure (ESM Atlas)",
        command=lambda: run_job("3D Structure Prediction", lambda: pred_ops.predict_structure(raw_input())),
    ).pack(anchor="w", padx=5, pady=2)
    tk.Button(
        buttons,
        text="Find IDs and download matches",
        command=lambda: run_job(
            "ID lookup + download",
            lambda: pred_ops.find_and_download_matches(raw_input()),
        ),
    ).pack(anchor="w", padx=5, pady=2)
    tk.Button(
        buttons,
        text="Fetch AlphaFold DB (sequence or UniProt ID)",
        command=lambda: run_job("AlphaFold DB", lambda: pred_ops.fetch_alphafold_db(raw_input())),
    ).pack(anchor="w", padx=5, pady=2)
    tk.Button(
        buttons,
        text="Fetch RCSB structure (sequence or PDB ID)",
        command=lambda: run_job("RCSB PDB", lambda: pred_ops.fetch_rcsb_pdb(raw_input())),
    ).pack(anchor="w", padx=5, pady=2)
    tk.Button(buttons, text="Load Local PDB", command=load_pdb).pack(anchor="w", padx=5, pady=2)
    tk.Button(
        buttons,
        text="Open in BioViewer",
        command=lambda: show("Open PDB", pred_ops.open_in_bioviewer(pred_ops.DEFAULT_PDB_NAME)),
    ).pack(anchor="w", padx=5, pady=2)
    tk.Button(
        buttons,
        text="Predict Secondary Structure",
        command=lambda: show("Secondary Structure", pred_ops.format_ss_report(raw_input())),
    ).pack(anchor="w", padx=5, pady=2)
    tk.Button(
        buttons,
        text="Calculate Hydrophobicity",
        command=lambda: show("Hydrophobicity", pred_ops.format_hydrophobicity_report(raw_input())),
    ).pack(anchor="w", padx=5, pady=2)
    tk.Button(
        buttons,
        text="Calculate Molecular Weight",
        command=lambda: show("Molecular Weight", pred_ops.format_mw_report(raw_input())),
    ).pack(anchor="w", padx=5, pady=2)


if __name__ == "__main__":
    create_main_window()
