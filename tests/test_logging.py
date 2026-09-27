"""Logging configuration and the privacy rule for sequence data."""

from __future__ import annotations

import json
import logging

from bioseqinsight.services.logging_setup import (
    configure_logging,
    get_logger,
    log_operation,
    sequence_reference,
)

UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)


class TestSequenceReference:
    def test_reference_does_not_contain_the_sequence(self):
        reference = sequence_reference(UBIQUITIN)
        assert UBIQUITIN not in reference
        assert UBIQUITIN[:20] not in reference

    def test_reference_records_the_length(self):
        assert "len=76" in sequence_reference(UBIQUITIN)

    def test_reference_is_stable(self):
        assert sequence_reference(UBIQUITIN) == sequence_reference(UBIQUITIN)

    def test_different_sequences_differ(self):
        assert sequence_reference("AAAA") != sequence_reference("CCCC")

    def test_empty_input_is_safe(self):
        assert sequence_reference("").startswith("len=0")


class TestConfiguration:
    def test_file_handler_writes(self, tmp_path):
        log_file = tmp_path / "logs" / "test.log"
        logger = configure_logging(log_file=log_file, console=False)
        logger.info("hello")
        for handler in logger.handlers:
            handler.flush()
        assert log_file.is_file() and "hello" in log_file.read_text()

    def test_json_format_emits_one_object_per_line(self, tmp_path):
        log_file = tmp_path / "test.log"
        logger = configure_logging(log_file=log_file, console=False, json_format=True)
        log_operation(logger, "retrieve", resource="rcsb", status="ok", elapsed_s=1.234)
        for handler in logger.handlers:
            handler.flush()
        record = json.loads(log_file.read_text().strip().splitlines()[-1])
        assert record["operation"] == "retrieve"
        assert record["resource"] == "rcsb"
        assert record["elapsed_s"] == 1.234

    def test_reconfiguring_does_not_duplicate_handlers(self, tmp_path):
        configure_logging(log_file=tmp_path / "a.log", console=False)
        logger = configure_logging(log_file=tmp_path / "a.log", console=False)
        assert len(logger.handlers) == 1

    def test_child_loggers_are_namespaced(self):
        assert get_logger("rcsb").name == "bioseqinsight.rcsb"

    def test_level_is_applied(self, tmp_path):
        logger = configure_logging(log_file=tmp_path / "a.log", level="WARNING", console=False)
        assert logger.level == logging.WARNING
