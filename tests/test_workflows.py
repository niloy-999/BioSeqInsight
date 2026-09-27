"""Batch analysis, project management and exporters."""

from __future__ import annotations

import csv
import io
import json
import zipfile
from pathlib import Path

import pytest

from bioseqinsight.core.alphabet import FastaRecord
from bioseqinsight.io.export import (
    BATCH_COLUMNS,
    build_html_report,
    rows_to_delimited,
    write_csv,
    write_html_report,
    write_json,
    write_tsv,
)
from bioseqinsight.models.results import MappingLevel
from bioseqinsight.services.cache import ResponseCache
from bioseqinsight.structures.manager import StructureManager
from bioseqinsight.workflows.batch import BatchOptions, BatchRunner
from bioseqinsight.workflows.project import Project, ProjectError, list_projects
from conftest import UBIQUITIN, make_client

DNA = "ATG" + "GCTAGCTTACCGATGGC" * 4 + "TAA"


def records() -> list[FastaRecord]:
    return [
        FastaRecord("dna1", "test nucleotide", DNA),
        FastaRecord("prot1", "ubiquitin", UBIQUITIN),
    ]


class TestBatchLocal:
    def test_both_alphabets_are_detected(self, tmp_settings):
        rows, summary = BatchRunner(tmp_settings).run(records())
        assert [row.input_type for row in rows] == ["dna", "protein"]
        assert summary.total == 2 and summary.failed == 0

    def test_dna_row_carries_a_sequence_result(self, tmp_settings):
        rows, _ = BatchRunner(tmp_settings).run(records())
        assert rows[0].sequence_result is not None
        assert rows[0].protein_result is None

    def test_protein_row_carries_a_protein_result(self, tmp_settings):
        rows, _ = BatchRunner(tmp_settings).run(records())
        assert rows[1].protein_result is not None
        assert rows[1].protein_result.length == 76

    def test_a_bad_record_does_not_abort_the_run(self, tmp_settings):
        bad = records() + [FastaRecord("broken", "", "@@@@@@@@@@@@")]
        rows, summary = BatchRunner(tmp_settings).run(bad)
        assert summary.total == 3
        assert summary.failed == 1
        assert rows[0].sequence_result is not None  # the good rows still ran

    def test_forced_alphabet_skips_detection(self, tmp_settings):
        rows, _ = BatchRunner(tmp_settings).run(
            [FastaRecord("x", "", "ACGTACGTACGTACGT")],
            BatchOptions(force_alphabet="protein"),
        )
        assert rows[0].input_type == "protein"

    def test_progress_callback_is_called_for_every_record(self, tmp_settings):
        seen = []
        BatchRunner(tmp_settings).run(
            records(), progress=lambda done, total, message: seen.append((done, total))
        )
        assert seen == [(1, 2), (2, 2)]

    def test_empty_input_is_handled(self, tmp_settings):
        rows, summary = BatchRunner(tmp_settings).run([])
        assert rows == [] and summary.total == 0

    def test_summary_reports_throughput(self, tmp_settings):
        _, summary = BatchRunner(tmp_settings).run(records())
        assert summary.to_dict()["sequences_per_second"] > 0


class TestBatchWithStructures:
    def test_structures_are_attached_and_validated(self, tmp_settings, working_transport):
        manager = StructureManager(
            tmp_settings, make_client(working_transport), ResponseCache(tmp_settings.cache_dir)
        )
        runner = BatchRunner(tmp_settings, manager=manager)
        rows, summary = runner.run(
            [FastaRecord("ubq", "", UBIQUITIN)],
            BatchOptions(include_structures=True, workers=1),
        )
        assert summary.structures_retrieved == 1
        assert rows[0].structures
        assert rows[0].best_mapping in {"M3", "M4"}

    def test_offline_batch_still_produces_local_results(self, tmp_settings):
        tmp_settings.offline = True
        rows, summary = BatchRunner(tmp_settings).run(
            records(), BatchOptions(include_structures=True, workers=1)
        )
        assert rows[1].protein_result is not None
        assert summary.structures_retrieved == 0

    def test_mapping_counts_are_tallied(self, tmp_settings, working_transport):
        manager = StructureManager(
            tmp_settings, make_client(working_transport), ResponseCache(tmp_settings.cache_dir)
        )
        rows, summary = BatchRunner(tmp_settings, manager=manager).run(
            [FastaRecord("ubq", "", UBIQUITIN)],
            BatchOptions(include_structures=True, workers=1),
        )
        assert sum(summary.mapping_counts.values()) == 1


class TestExport:
    def test_csv_has_the_documented_columns(self, tmp_settings):
        rows, _ = BatchRunner(tmp_settings).run(records())
        reader = csv.DictReader(io.StringIO(rows_to_delimited(rows)))
        assert reader.fieldnames == list(BATCH_COLUMNS)

    def test_csv_has_one_row_per_record(self, tmp_settings):
        rows, _ = BatchRunner(tmp_settings).run(records())
        parsed = list(csv.DictReader(io.StringIO(rows_to_delimited(rows))))
        assert len(parsed) == 2
        assert parsed[0]["identifier"] == "dna1"

    def test_tsv_uses_tabs(self, tmp_settings, tmp_path):
        rows, _ = BatchRunner(tmp_settings).run(records())
        path = write_tsv(rows, tmp_path / "out.tsv")
        assert "\t" in Path(path).read_text().splitlines()[0]

    def test_rows_without_structures_are_reported_as_m0(self, tmp_settings):
        rows, _ = BatchRunner(tmp_settings).run(records())
        parsed = list(csv.DictReader(io.StringIO(rows_to_delimited(rows))))
        assert parsed[0]["mapping_level"] == MappingLevel.M0_NO_RESULT.value

    def test_json_export_round_trips(self, tmp_settings, tmp_path):
        rows, summary = BatchRunner(tmp_settings).run(records())
        path = write_json({"rows": [r.to_dict() for r in rows]}, tmp_path / "out.json")
        assert len(json.loads(Path(path).read_text())["rows"]) == 2

    def test_html_report_is_self_contained(self, tmp_settings):
        rows, _ = BatchRunner(tmp_settings).run(records())
        html = build_html_report(rows, title="Test report")
        assert html.startswith("<!DOCTYPE html>")
        assert "Test report" in html
        assert "<style>" in html  # CSS inlined, no external stylesheet
        assert "src=" not in html  # no external scripts or images

    def test_html_report_explains_the_mapping_levels(self, tmp_settings):
        rows, _ = BatchRunner(tmp_settings).run(records())
        html = build_html_report(rows)
        assert "M4" in html and "M0" in html

    def test_html_escapes_user_supplied_identifiers(self, tmp_settings):
        rows, _ = BatchRunner(tmp_settings).run(
            [FastaRecord("<script>alert(1)</script>", "", DNA)]
        )
        assert "<script>alert(1)</script>" not in build_html_report(rows)

    def test_csv_file_is_written(self, tmp_settings, tmp_path):
        rows, _ = BatchRunner(tmp_settings).run(records())
        assert Path(write_csv(rows, tmp_path / "a" / "out.csv")).is_file()

    def test_html_file_is_written(self, tmp_settings, tmp_path):
        rows, _ = BatchRunner(tmp_settings).run(records())
        assert Path(write_html_report(rows, tmp_path / "out.html")).is_file()


class TestProject:
    def test_create_makes_the_expected_layout(self, tmp_path, tmp_settings):
        project = Project.create(tmp_path / "p", name="demo", settings=tmp_settings)
        for name in ("sequences", "structures", "results", "logs"):
            assert (project.path / name).is_dir()
        assert (project.path / "project.json").is_file()
        assert (project.path / "metadata.json").is_file()

    def test_settings_snapshot_is_stored(self, tmp_path, tmp_settings):
        project = Project.create(tmp_path / "p", settings=tmp_settings)
        metadata = json.loads((project.path / "metadata.json").read_text())
        assert metadata["settings"]["max_attempts"] == tmp_settings.max_attempts
        assert metadata["software_version"]

    def test_creating_over_a_non_empty_directory_is_refused(self, tmp_path, tmp_settings):
        target = tmp_path / "p"
        target.mkdir()
        (target / "something.txt").write_text("x")
        with pytest.raises(ProjectError, match="not empty"):
            Project.create(target, settings=tmp_settings)

    def test_opening_a_non_project_is_refused(self, tmp_path, tmp_settings):
        with pytest.raises(ProjectError, match="not a BioSeqInsight project"):
            Project.open(tmp_path)

    def test_future_format_version_is_refused(self, tmp_path, tmp_settings):
        project = Project.create(tmp_path / "p", settings=tmp_settings)
        manifest = json.loads((project.path / "project.json").read_text())
        manifest["format_version"] = 99
        (project.path / "project.json").write_text(json.dumps(manifest))
        with pytest.raises(ProjectError, match="format version"):
            Project.open(project.path)

    def test_sequences_round_trip(self, tmp_path, tmp_settings):
        project = Project.create(tmp_path / "p", settings=tmp_settings)
        project.add_sequences(records())
        reopened = Project.open(project.path)
        assert len(reopened.read_sequences()) == 2
        assert reopened.manifest.record_count == 2

    def test_results_are_saved_in_every_requested_format(self, tmp_path, tmp_settings):
        project = Project.create(tmp_path / "p", settings=tmp_settings)
        rows, _ = BatchRunner(tmp_settings).run(records())
        written = project.save_results(rows, formats=("csv", "json", "html"))
        assert set(written) == {"csv", "json", "html"}
        assert all(Path(path).is_file() for path in written.values())

    def test_export_and_import_round_trip(self, tmp_path, tmp_settings):
        project = Project.create(tmp_path / "p", name="demo", settings=tmp_settings)
        project.add_sequences(records())
        archive = project.export_zip(tmp_path / "demo.zip")
        assert zipfile.is_zipfile(archive)

        restored = Project.import_zip(archive, tmp_path / "restored")
        assert restored.manifest.name == "demo"
        assert len(restored.read_sequences()) == 2

    def test_import_refuses_unsafe_paths(self, tmp_path, tmp_settings):
        archive = tmp_path / "evil.zip"
        with zipfile.ZipFile(archive, "w") as zf:
            zf.writestr("../escape.txt", "x")
        with pytest.raises(ProjectError, match="unsafe path"):
            Project.import_zip(archive, tmp_path / "out")

    def test_notes_are_timestamped(self, tmp_path, tmp_settings):
        project = Project.create(tmp_path / "p", settings=tmp_settings)
        project.write_note("benchmark run 1")
        assert Project.open(project.path).manifest.notes[0].endswith("benchmark run 1")

    def test_listing_projects(self, tmp_path, tmp_settings):
        for name in ("a", "b"):
            Project.create(tmp_path / name, settings=tmp_settings)
        assert len(list_projects(tmp_path)) == 2

    def test_listing_a_missing_root_is_empty(self, tmp_path):
        assert list_projects(tmp_path / "absent") == []

    def test_import_fasta_from_disk(self, tmp_path, tmp_settings):
        source = tmp_path / "in.fasta"
        source.write_text(">x\nATGCATGC\n")
        project = Project.create(tmp_path / "p", settings=tmp_settings)
        project.import_fasta(source)
        assert project.read_sequences()[0].sequence == "ATGCATGC"

    def test_import_missing_fasta_is_reported(self, tmp_path, tmp_settings):
        project = Project.create(tmp_path / "p", settings=tmp_settings)
        with pytest.raises(ProjectError, match="No such FASTA"):
            project.import_fasta(tmp_path / "absent.fasta")
