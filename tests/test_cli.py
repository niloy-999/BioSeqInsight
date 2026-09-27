"""Command-line interface.

The CLI is how the manuscript's evaluation is reproduced, so it is tested
like any other interface: exit codes, output shape and error handling.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from bioseqinsight.cli import main

UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)
DNA = "ATG" + "GCTAGCTTACCGATGGC" * 4 + "TAA"


def run(argv, tmp_path):
    """Run the CLI with every path confined to a temporary directory."""
    return main(
        [
            "--quiet",
            "--offline",
            "--download-dir",
            str(tmp_path / "structures"),
            *argv,
        ]
    )


class TestDNACommand:
    def test_exit_code_is_zero(self, tmp_path):
        assert run(["dna", DNA], tmp_path) == 0

    def test_json_output_is_parseable(self, tmp_path):
        out = tmp_path / "r.json"
        run(["dna", DNA, "--json", "--out", str(out)], tmp_path)
        payload = json.loads(out.read_text())
        assert payload["length"] == len(DNA)
        assert payload["software_version"]

    def test_text_output_contains_the_key_numbers(self, tmp_path):
        out = tmp_path / "r.txt"
        run(["dna", DNA, "--out", str(out)], tmp_path)
        text = out.read_text()
        assert "GC content" in text and "Tm" in text

    def test_motif_option(self, tmp_path):
        out = tmp_path / "r.txt"
        run(["dna", DNA, "--motif", "GCTAGC", "--out", str(out)], tmp_path)
        assert "Motif GCTAGC" in out.read_text()

    def test_reads_a_fasta_file(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(f">seq1\n{DNA}\n")
        out = tmp_path / "r.json"
        run(["dna", "--file", str(fasta), "--json", "--out", str(out)], tmp_path)
        assert json.loads(out.read_text())["identifier"] == "seq1"

    def test_missing_file_is_a_user_error(self, tmp_path):
        assert run(["dna", "--file", str(tmp_path / "absent.fasta")], tmp_path) == 2

    def test_invalid_sequence_is_a_user_error(self, tmp_path):
        assert run(["dna", "@@@@@@"], tmp_path) == 2

    def test_genetic_code_table_option(self, tmp_path):
        out = tmp_path / "r.json"
        assert run(["dna", DNA, "--table", "11", "--json", "--out", str(out)], tmp_path) == 0


class TestProteinCommand:
    def test_exit_code_and_output(self, tmp_path):
        out = tmp_path / "p.txt"
        assert run(["protein", UBIQUITIN, "--out", str(out)], tmp_path) == 0
        assert "Molecular weight" in out.read_text()

    def test_json_output(self, tmp_path):
        out = tmp_path / "p.json"
        run(["protein", UBIQUITIN, "--json", "--out", str(out)], tmp_path)
        assert json.loads(out.read_text())["length"] == 76

    def test_sketch_flag(self, tmp_path):
        out = tmp_path / "p.txt"
        run(["protein", UBIQUITIN, "--sketch", "--out", str(out)], tmp_path)
        assert "propensity" in out.read_text().lower()

    def test_empty_input_is_a_user_error(self, tmp_path):
        assert run(["protein", "   "], tmp_path) == 2


class TestBatchCommand:
    def test_batch_writes_csv_and_json(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(f">dna1\n{DNA}\n>prot1\n{UBIQUITIN}\n")
        prefix = tmp_path / "out"
        assert run(
            ["batch", str(fasta), "--out-prefix", str(prefix), "--formats", "csv,json"], tmp_path
        ) == 0
        assert Path(f"{prefix}.csv").is_file()
        payload = json.loads(Path(f"{prefix}.json").read_text())
        assert payload["summary"]["total"] == 2

    def test_batch_writes_html(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(f">prot1\n{UBIQUITIN}\n")
        prefix = tmp_path / "out"
        run(["batch", str(fasta), "--out-prefix", str(prefix), "--formats", "html"], tmp_path)
        assert "<!DOCTYPE html>" in Path(f"{prefix}.html").read_text()

    def test_missing_input_is_a_user_error(self, tmp_path):
        assert run(["batch", str(tmp_path / "absent.fasta")], tmp_path) == 2

    def test_empty_fasta_is_a_user_error(self, tmp_path):
        empty = tmp_path / "empty.fasta"
        empty.write_text("")
        assert run(["batch", str(empty)], tmp_path) == 2

    def test_batch_into_a_project(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(f">prot1\n{UBIQUITIN}\n")
        project = tmp_path / "proj"
        assert run(["batch", str(fasta), "--project", str(project)], tmp_path) == 0
        assert (project / "project.json").is_file()
        assert list((project / "results").glob("*.csv"))


class TestProjectCommand:
    def test_create_info_and_export(self, tmp_path):
        project = tmp_path / "p"
        assert run(["project", "create", str(project), "--name", "demo"], tmp_path) == 0
        assert run(["project", "info", str(project)], tmp_path) == 0
        assert run(["project", "export", str(project), "--out", str(tmp_path / "p.zip")], tmp_path) == 0
        assert (tmp_path / "p.zip").is_file()

    def test_import_round_trip(self, tmp_path):
        project = tmp_path / "p"
        run(["project", "create", str(project), "--name", "demo"], tmp_path)
        run(["project", "export", str(project), "--out", str(tmp_path / "p.zip")], tmp_path)
        assert run(["project", "import", str(tmp_path / "p.zip"), str(tmp_path / "restored")], tmp_path) == 0

    def test_info_on_a_non_project_is_a_user_error(self, tmp_path):
        assert run(["project", "info", str(tmp_path)], tmp_path) == 2

    def test_create_with_fasta(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(f">x\n{DNA}\n")
        assert run(["project", "create", str(tmp_path / "p"), "--fasta", str(fasta)], tmp_path) == 0


class TestOtherCommands:
    def test_version_reports_the_environment(self, tmp_path):
        assert run(["version"], tmp_path) == 0

    def test_config_prints_settings(self, tmp_path):
        assert run(["config"], tmp_path) == 0

    def test_config_save(self, tmp_path):
        target = tmp_path / "settings.json"
        assert run(["config", "--save", str(target)], tmp_path) == 0
        assert json.loads(target.read_text())["max_attempts"] >= 1

    def test_cache_stats(self, tmp_path):
        assert run(["cache", "stats"], tmp_path) == 0

    def test_cache_clear(self, tmp_path):
        assert run(["cache", "clear"], tmp_path) == 0

    def test_structure_in_offline_mode_reports_no_result(self, tmp_path):
        # Offline, every provider is skipped, so the command reports failure
        # rather than pretending a structure was found.
        assert run(["structure", "P24941"], tmp_path) == 3

    def test_unknown_command_exits_with_usage_error(self, tmp_path):
        with pytest.raises(SystemExit):
            main(["not-a-command"])
