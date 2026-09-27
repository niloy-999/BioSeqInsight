#!/usr/bin/env python3
"""Integrated versus manual workflow comparison.

The editorial criticism that v1.0 offered no demonstrated "added value over
existing tools" needs an answer that does not rest on assertion. This script
provides the part of that answer which can be derived objectively.

What is measured here
---------------------
**Interaction cost.** For each of five representative tasks, the manual
protocol (documented in ``benchmarks/protocols/manual_workflow_protocol.md``)
is expressed as an explicit list of steps. Each step is labelled with the
application it happens in and whether it involves transcribing data by hand
between tools. Counting steps, distinct applications and manual data
transfers is objective: two people following the protocol will count the
same, and the protocol is published so the count can be checked.

**Machine time.** For the same tasks, the wall-clock time BioSeqInsight takes
to produce the equivalent output is measured directly, offline.

What is *not* measured here
---------------------------
Human completion time, error rates and subjective workload. Those require
participants, and inventing them would be worse than omitting them. They are
the subject of ``benchmarks/protocols/usability_study_protocol.md``, which is
designed to be run by the authors and reported alongside these counts.

The distinction matters: interaction counts and machine timings are facts
about the software; task times and error rates are facts about people using
it. This script reports only the former.

Usage
-----
    python benchmarks/scripts/run_workflow_comparison.py
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from bioseqinsight import __version__  # noqa: E402
from bioseqinsight.config.settings import Settings  # noqa: E402
from bioseqinsight.core.alphabet import FastaRecord  # noqa: E402
from bioseqinsight.core.protein import ProteinAnalyzer  # noqa: E402
from bioseqinsight.core.sequence import SequenceAnalyzer  # noqa: E402
from bioseqinsight.io.export import write_csv  # noqa: E402
from bioseqinsight.workflows.batch import BatchOptions, BatchRunner  # noqa: E402

UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)
DEMO_DNA = "ATG" + "GCTAGCTTACCGATGGCATTACGGATCCTAGGCTTAACCGGATCAGCT" * 6 + "TAA"


@dataclass
class Step:
    """One action a user performs, in one application."""

    description: str
    application: str
    manual_transfer: bool = False  # copying or retyping data between tools


@dataclass
class Task:
    name: str
    goal: str
    manual: list[Step]
    integrated: list[Step]
    machine_time_s: float | None = None
    notes: str = ""
    extra: dict = field(default_factory=dict)

    def counts(self, steps: list[Step]) -> dict:
        return {
            "steps": len(steps),
            "applications": len(sorted({step.application for step in steps})),
            "application_list": sorted({step.application for step in steps}),
            "manual_data_transfers": sum(1 for step in steps if step.manual_transfer),
        }

    def to_dict(self) -> dict:
        manual = self.counts(self.manual)
        integrated = self.counts(self.integrated)
        # Where a task's manual protocol contains a per-sequence loop, the
        # expanded count is the honest one to compare against.
        manual_steps = self.extra.get("manual_expanded_steps", manual["steps"])
        return {
            "task": self.name,
            "goal": self.goal,
            "manual": manual,
            "integrated": integrated,
            "step_reduction_percent": _reduction(manual_steps, integrated["steps"]),
            "application_reduction": manual["applications"] - integrated["applications"],
            "manual_transfers_eliminated": (
                manual["manual_data_transfers"] - integrated["manual_data_transfers"]
            ),
            "bioseqinsight_machine_time_s": (
                round(self.machine_time_s, 3) if self.machine_time_s is not None else None
            ),
            "notes": self.notes,
            **self.extra,
        }


def _reduction(before: int, after: int) -> float:
    return round(100.0 * (before - after) / before, 1) if before else 0.0


# ---------------------------------------------------------------------------
# Task definitions. The manual protocols mirror
# benchmarks/protocols/manual_workflow_protocol.md exactly.
# ---------------------------------------------------------------------------


def build_tasks() -> list[Task]:
    return [
        Task(
            name="T1_single_sequence_characterisation",
            goal=(
                "For one coding sequence: GC content, melting temperature, longest ORF "
                "and the physicochemical properties of its product."
            ),
            manual=[
                Step("Open a GC-content web tool", "web tool A"),
                Step("Paste the sequence", "web tool A", manual_transfer=True),
                Step("Read and record GC content", "notes", manual_transfer=True),
                Step("Open a melting-temperature calculator", "web tool B"),
                Step("Paste the sequence again", "web tool B", manual_transfer=True),
                Step("Set primer and salt concentrations", "web tool B"),
                Step("Read and record Tm", "notes", manual_transfer=True),
                Step("Open an ORF finder", "web tool C"),
                Step("Paste the sequence again", "web tool C", manual_transfer=True),
                Step("Choose the genetic code and minimum length", "web tool C"),
                Step("Read the ORF table", "web tool C"),
                Step("Copy the longest ORF peptide", "web tool C", manual_transfer=True),
                Step("Open a protein parameter tool", "web tool D"),
                Step("Paste the peptide", "web tool D", manual_transfer=True),
                Step("Read molecular weight, pI and GRAVY", "web tool D"),
                Step("Transcribe all values into a notebook or spreadsheet", "notes", manual_transfer=True),
            ],
            integrated=[
                Step("Paste the sequence into the Sequence tab", "BioSeqInsight", manual_transfer=True),
                Step("Click Full analysis", "BioSeqInsight"),
                Step("Click Export CSV", "BioSeqInsight"),
            ],
            notes=(
                "The manual route recomputes nothing: it retypes the same sequence into "
                "four tools, each with its own defaults, and the analysis conditions are "
                "not recorded anywhere."
            ),
        ),
        Task(
            name="T2_structure_retrieval_with_identity_check",
            goal=(
                "Obtain a structure for a protein and establish whether it really is "
                "that protein rather than a homologue or a different molecule."
            ),
            manual=[
                Step("Open the UniProt entry", "UniProt"),
                Step("Copy the canonical sequence", "UniProt", manual_transfer=True),
                Step("Open the RCSB search page", "RCSB PDB"),
                Step("Paste the sequence and run a sequence search", "RCSB PDB", manual_transfer=True),
                Step("Inspect the hit list and choose an entry", "RCSB PDB"),
                Step("Open the entry page and check its UniProt cross-reference", "RCSB PDB"),
                Step("Download the coordinate file", "RCSB PDB"),
                Step("Extract the chain sequence from the file", "text editor"),
                Step("Open a pairwise alignment tool", "EMBOSS/BLAST"),
                Step("Paste the query sequence", "EMBOSS/BLAST", manual_transfer=True),
                Step("Paste the chain sequence", "EMBOSS/BLAST", manual_transfer=True),
                Step("Run the alignment and read identity and coverage", "EMBOSS/BLAST"),
                Step("Decide whether the match is acceptable", "judgement"),
                Step("If not, open the AlphaFold DB entry", "AlphaFold DB"),
                Step("Download the predicted model", "AlphaFold DB"),
                Step("Read the mean pLDDT from the file or the page", "AlphaFold DB"),
                Step("Record the decision and the evidence for it", "notes", manual_transfer=True),
            ],
            integrated=[
                Step("Enter the accession in the Structure tab", "BioSeqInsight", manual_transfer=True),
                Step("Click Retrieve and validate", "BioSeqInsight"),
                Step("Read the mapping level, identity and coverage", "BioSeqInsight"),
            ],
            notes=(
                "The manual alignment step is the one most often skipped, and skipping "
                "it is exactly how v1.0's benchmark reported 25/25 successful "
                "retrievals while 7 of them were the wrong molecule."
            ),
        ),
        Task(
            name="T3_batch_of_fifty_sequences",
            goal="Produce one table of properties for fifty sequences.",
            manual=[
                Step("Split the FASTA file into individual sequences", "text editor"),
                Step("For each sequence, repeat task T1 (16 steps)", "web tools A-D", manual_transfer=True),
                Step("Paste each result into a spreadsheet", "spreadsheet", manual_transfer=True),
                Step("Check for transcription errors", "spreadsheet"),
                Step("Format the table", "spreadsheet"),
            ],
            integrated=[
                Step("Choose the FASTA file in the Batch tab", "BioSeqInsight"),
                Step("Click Run batch", "BioSeqInsight"),
                Step("Click Export CSV", "BioSeqInsight"),
            ],
            notes=(
                "The manual step count below expands the per-sequence loop: 50 x 16 "
                "steps plus the surrounding spreadsheet work. The integrated count "
                "does not grow with the number of sequences, which is the point."
            ),
            extra={"manual_expanded_steps": 50 * 16 + 4},
        ),
        Task(
            name="T4_reproduce_an_earlier_analysis",
            goal="Re-run an analysis from six months ago and get the same numbers.",
            manual=[
                Step("Find the notebook entry or spreadsheet", "notes"),
                Step("Identify which web tools were used", "notes"),
                Step("Check whether each tool still exists and behaves the same", "web tools"),
                Step("Guess the parameters that were used", "judgement"),
                Step("Repeat the whole analysis", "web tools A-D", manual_transfer=True),
                Step("Compare with the recorded values and reconcile differences", "spreadsheet"),
            ],
            integrated=[
                Step("Open the project folder", "BioSeqInsight"),
                Step("Read the stored settings snapshot and software version", "BioSeqInsight"),
                Step("Re-run the batch", "BioSeqInsight"),
            ],
            notes=(
                "Each BioSeqInsight project stores the software version, the complete "
                "settings and the input sequences, so the run is re-executable. The "
                "manual route depends on what the analyst happened to write down."
            ),
        ),
        Task(
            name="T5_share_an_analysis_with_a_collaborator",
            goal="Give a colleague everything needed to check and continue the work.",
            manual=[
                Step("Collect the input files", "file manager"),
                Step("Collect the downloaded structures", "file manager"),
                Step("Collect the spreadsheet of results", "file manager"),
                Step("Write down which tools and settings were used", "notes", manual_transfer=True),
                Step("Zip the folder by hand", "file manager"),
                Step("Send it and explain the layout", "email"),
            ],
            integrated=[
                Step("Project > Export project as zip", "BioSeqInsight"),
                Step("Send the archive", "email"),
            ],
            notes=(
                "The exported archive contains the sequences, the retrieved "
                "coordinates, the results, the logs, the settings snapshot and the "
                "software version, and re-imports into a working project."
            ),
        ),
    ]


# ---------------------------------------------------------------------------
# Machine timings for the integrated route
# ---------------------------------------------------------------------------


def measure_machine_times(tasks: list[Task]) -> None:
    by_name = {task.name: task for task in tasks}
    with TemporaryDirectory() as tmp:
        settings = Settings.from_dict(
            {
                "download_dir": f"{tmp}/s",
                "cache_dir": f"{tmp}/c",
                "projects_dir": f"{tmp}/p",
                "log_file": f"{tmp}/l/x.log",
                "offline": True,
            }
        )

        started = time.perf_counter()
        SequenceAnalyzer(DEMO_DNA).analyze(min_orf_aa=30)
        ProteinAnalyzer(UBIQUITIN).analyze()
        by_name["T1_single_sequence_characterisation"].machine_time_s = (
            time.perf_counter() - started
        )

        records = [FastaRecord(f"seq{i}", "", DEMO_DNA) for i in range(50)]
        runner = BatchRunner(settings)
        started = time.perf_counter()
        rows, _ = runner.run(records, BatchOptions(include_structures=False))
        write_csv(rows, f"{tmp}/batch.csv")
        by_name["T3_batch_of_fifty_sequences"].machine_time_s = time.perf_counter() - started

        from bioseqinsight.workflows.project import Project

        project = Project.create(f"{tmp}/project", name="comparison", settings=settings)
        project.add_sequences(records)
        project.save_results(rows, formats=("csv",))
        started = time.perf_counter()
        project.export_zip(f"{tmp}/project.zip")
        by_name["T5_share_an_analysis_with_a_collaborator"].machine_time_s = (
            time.perf_counter() - started
        )

        started = time.perf_counter()
        reopened = Project.open(f"{tmp}/project")
        runner.run(reopened.read_sequences(), BatchOptions(include_structures=False))
        by_name["T4_reproduce_an_earlier_analysis"].machine_time_s = time.perf_counter() - started

    # T2 depends on external services; its machine time is reported by the
    # live structure benchmark, not invented here.
    by_name["T2_structure_retrieval_with_identity_check"].notes += (
        " Machine time for this task is network-bound and is reported by "
        "run_structure_benchmark.py."
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="benchmarks/results/workflow_comparison.json")
    args = parser.parse_args(argv)

    tasks = build_tasks()
    measure_machine_times(tasks)
    entries = [task.to_dict() for task in tasks]

    total_manual = sum(
        entry.get("manual_expanded_steps", entry["manual"]["steps"]) for entry in entries
    )
    total_integrated = sum(entry["integrated"]["steps"] for entry in entries)

    report = {
        "software_version": __version__,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "method": (
            "Step, application and manual-data-transfer counts derived from the "
            "published manual protocol (benchmarks/protocols/manual_workflow_protocol.md). "
            "Machine timings measured offline on this host. Human completion times and "
            "error rates are NOT included here; see usability_study_protocol.md."
        ),
        "tasks": entries,
        "totals": {
            "manual_steps": total_manual,
            "integrated_steps": total_integrated,
            "step_reduction_percent": _reduction(total_manual, total_integrated),
            "distinct_applications_manual": len(
                sorted({app for task in tasks for step in task.manual for app in [step.application]})
            ),
            "distinct_applications_integrated": len(
                sorted(
                    {app for task in tasks for step in task.integrated for app in [step.application]}
                )
            ),
            "manual_data_transfers_manual": sum(
                entry["manual"]["manual_data_transfers"] for entry in entries
            ),
            "manual_data_transfers_integrated": sum(
                entry["integrated"]["manual_data_transfers"] for entry in entries
            ),
        },
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")

    width = 44
    print(f"{'task':<{width}}{'manual':>8}{'integrated':>12}{'reduction':>11}{'machine s':>11}")
    print("-" * (width + 42))
    for entry in entries:
        manual_steps = entry.get("manual_expanded_steps", entry["manual"]["steps"])
        machine = entry["bioseqinsight_machine_time_s"]
        print(
            f"{entry['task']:<{width}}{manual_steps:>8}{entry['integrated']['steps']:>12}"
            f"{entry['step_reduction_percent']:>10.0f}%"
            f"{(f'{machine:.3f}' if machine is not None else 'network'):>11}"
        )
    print("-" * (width + 42))
    print(json.dumps(report["totals"], indent=2))
    print(f"\nWritten to {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
