"""The partial-structure qualifier (Amendment 7).

A structure that IS the queried protein but covers only part of it (a domain
construct) must be reported as such, without moving the mapping level and
without ever becoming M4.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from bioseqinsight.config.settings import Settings
from bioseqinsight.models.results import MappingLevel
from bioseqinsight.structures.pdbio import parse_pdb
from bioseqinsight.structures.validation import (
    classify,
    describe,
    is_partial_structure,
    validate_structure,
)

DATA = Path(__file__).parent / "data" / "1ubq_excerpt.pdb"


@pytest.fixture()
def settings() -> Settings:
    return Settings()


@pytest.mark.parametrize(
    "identity, coverage, accession_match, expected",
    [
        (99.9, 27.0, True, True),  # the B2 fragment pattern
        (99.9, 79.9, True, True),  # just below the high-coverage threshold
        (99.9, 80.0, True, False),  # at the threshold: not partial
        (99.9, 27.0, False, False),  # accession not confirmed
        (99.9, 27.0, None, False),  # no accession queried
        (94.9, 27.0, True, False),  # identity below the high threshold
        (None, 27.0, True, False),  # no comparison was possible
        (99.9, None, True, False),
    ],
)
def test_is_partial_structure(settings, identity, coverage, accession_match, expected):
    assert is_partial_structure(identity, coverage, accession_match, settings) is expected


def test_partial_flag_does_not_change_the_mapping_level(settings):
    # The same inputs that are flagged partial are still M1 (coverage < 40):
    # the qualifier is additive and the scale is untouched.
    assert classify(99.9, 27.0, True, settings) is MappingLevel.M1_UNVERIFIED
    assert classify(99.9, 60.0, True, settings) is MappingLevel.M2_SEQUENCE_MATCH


def test_partial_message_is_accurate_and_not_a_mismatch_warning():
    text = describe(
        MappingLevel.M1_UNVERIFIED, 99.9, 27.0, True, "P0CG48", ["P0CG48"], partial=True
    )
    assert text.startswith("Partial structure of the requested protein")
    assert "P0CG48" in text and "27.0%" in text
    assert "does not correspond" not in text


def test_non_partial_m1_still_warns():
    text = describe(MappingLevel.M1_UNVERIFIED, 20.0, 60.0, False, "P0CG48", ["Q99999"])
    assert text.startswith("WARNING")


def test_partial_never_affects_exact_matches():
    text = describe(MappingLevel.M4_EXACT, 100.0, 100.0, True, "P0CG48", ["P0CG48"], partial=True)
    assert text.startswith("Exact biological identity confirmed")


@pytest.mark.skipif(not DATA.is_file(), reason="structure excerpt not available")
class TestEndToEnd:
    """validate_structure on a real excerpt, with the query built from its own chain."""

    @staticmethod
    def _chain_sequence() -> str:
        structure = parse_pdb(DATA.read_text())
        assert structure.chains, "excerpt has no readable chain"
        return structure.chains[0].sequence

    def test_fragment_of_the_right_protein_is_flagged_partial(self, settings):
        chain = self._chain_sequence()
        query = chain + "W" * (3 * len(chain))  # the structure covers ~25% of this query
        report = validate_structure(
            DATA.read_text(),
            query_sequence=query,
            query_accession="P0CG48",
            settings=settings,
            extra_accessions=["P0CG48"],
        )
        assert report.coverage_percent is not None and report.coverage_percent < 80
        assert report.identity_percent is not None and report.identity_percent >= 95
        assert report.partial_structure is True
        assert report.mapping_level in (MappingLevel.M1_UNVERIFIED, MappingLevel.M2_SEQUENCE_MATCH)
        assert not report.exact_biological_identity
        assert report.message.startswith("Partial structure")

    def test_same_fragment_with_wrong_accession_is_not_partial(self, settings):
        chain = self._chain_sequence()
        report = validate_structure(
            DATA.read_text(),
            query_sequence=chain + "W" * (3 * len(chain)),
            query_accession="P99999",
            settings=settings,
        )
        assert report.partial_structure is False
        assert report.mapping_level is not MappingLevel.M4_EXACT

    def test_full_length_exact_match_is_not_partial(self, settings):
        chain = self._chain_sequence()
        report = validate_structure(
            DATA.read_text(),
            query_sequence=chain,
            query_accession="P0CG48",
            settings=settings,
            extra_accessions=["P0CG48"],
        )
        assert report.mapping_level is MappingLevel.M4_EXACT
        assert report.partial_structure is False

    def test_flag_is_exported(self, settings):
        chain = self._chain_sequence()
        report = validate_structure(
            DATA.read_text(),
            query_sequence=chain + "W" * (3 * len(chain)),
            query_accession="P0CG48",
            settings=settings,
            extra_accessions=["P0CG48"],
        )
        assert report.to_dict()["partial_structure"] is True
