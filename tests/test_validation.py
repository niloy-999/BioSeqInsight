"""The M0-M4 sequence-to-structure mapping taxonomy.

These tests encode the behaviour that distinguishes v2.0 from v1.0: a
retrieved structure is never reported as the queried molecule unless the
sequence and the accession both agree.
"""

from __future__ import annotations

import pytest

from bioseqinsight.models.results import MappingLevel
from bioseqinsight.structures.validation import (
    classify,
    describe,
    no_result_report,
    validate_structure,
)

UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)


class TestClassification:
    def test_exact_requires_accession_and_sequence(self, tmp_settings):
        assert classify(100.0, 100.0, True, tmp_settings) is MappingLevel.M4_EXACT

    def test_perfect_sequence_without_accession_is_only_high_confidence(self, tmp_settings):
        # This is the v1.0 failure mode: a perfect structural hit is not
        # proof that it is the requested biological entity.
        assert classify(100.0, 100.0, False, tmp_settings) is MappingLevel.M3_HIGH_CONFIDENCE

    def test_matching_accession_with_poor_sequence_is_not_exact(self, tmp_settings):
        assert classify(60.0, 60.0, True, tmp_settings) is MappingLevel.M2_SEQUENCE_MATCH

    def test_high_identity_low_coverage_is_not_exact(self, tmp_settings):
        assert classify(100.0, 50.0, True, tmp_settings) is MappingLevel.M2_SEQUENCE_MATCH

    def test_homologue_is_a_sequence_match(self, tmp_settings):
        assert classify(70.0, 90.0, False, tmp_settings) is MappingLevel.M2_SEQUENCE_MATCH

    def test_missing_measurements_give_unverified(self, tmp_settings):
        assert classify(None, None, True, tmp_settings) is MappingLevel.M1_UNVERIFIED

    def test_thresholds_are_configurable(self, tmp_settings):
        tmp_settings.exact_identity_threshold = 80.0
        tmp_settings.exact_coverage_threshold = 80.0
        assert classify(85.0, 85.0, True, tmp_settings) is MappingLevel.M4_EXACT

    def test_levels_are_ordered(self):
        ranks = [level.rank for level in MappingLevel]
        assert ranks == sorted(ranks)

    def test_every_level_has_a_label(self):
        assert all(level.label for level in MappingLevel)


class TestValidateStructure:
    def test_exact_match(self, ubiquitin_pdb, tmp_settings):
        report = validate_structure(
            ubiquitin_pdb,
            query_sequence=UBIQUITIN,
            query_accession="P62988",
            settings=tmp_settings,
        )
        assert report.mapping_level is MappingLevel.M4_EXACT
        assert report.identity_percent == pytest.approx(100.0)
        assert report.coverage_percent == pytest.approx(100.0)
        assert report.exact_biological_identity is True
        assert report.accession_match is True

    def test_wrong_accession_downgrades_a_perfect_sequence_hit(self, ubiquitin_pdb, tmp_settings):
        report = validate_structure(
            ubiquitin_pdb,
            query_sequence=UBIQUITIN,
            query_accession="P24941",  # CDK2, which 1UBQ is not
            settings=tmp_settings,
        )
        assert report.mapping_level is MappingLevel.M3_HIGH_CONFIDENCE
        assert report.accession_match is False
        assert "not confirmed" in report.message or "rather than" in report.message

    def test_unrelated_sequence_is_never_promoted_to_a_match(self, ubiquitin_pdb, tmp_settings):
        # A measured 0% identity must not be dressed up as a "related record".
        report = validate_structure(
            ubiquitin_pdb,
            query_sequence="W" * 76,
            query_accession=None,
            settings=tmp_settings,
        )
        assert report.mapping_level is MappingLevel.M1_UNVERIFIED
        assert report.identity_percent == pytest.approx(0.0)
        assert "WARNING" in report.message

    def test_homologue_lands_in_m2(self, ubiquitin_pdb, tmp_settings):
        # Mutate 20% of ubiquitin: still clearly related, not the same entry.
        mutated = "".join(
            ("A" if index % 5 == 0 else residue) for index, residue in enumerate(UBIQUITIN)
        )
        report = validate_structure(
            ubiquitin_pdb, query_sequence=mutated, query_accession=None, settings=tmp_settings
        )
        assert report.mapping_level is MappingLevel.M2_SEQUENCE_MATCH

    def test_missing_query_sequence_gives_unverified(self, ubiquitin_pdb, tmp_settings):
        report = validate_structure(
            ubiquitin_pdb, query_sequence=None, query_accession="P62988", settings=tmp_settings
        )
        assert report.mapping_level is MappingLevel.M1_UNVERIFIED
        assert "unverified" in report.message.lower()

    def test_payload_without_chains_gives_unverified(self, tmp_settings):
        report = validate_structure(
            "HEADER nothing\nEND\n",
            query_sequence=UBIQUITIN,
            query_accession="P62988",
            settings=tmp_settings,
        )
        assert report.mapping_level is MappingLevel.M1_UNVERIFIED

    def test_extra_accessions_from_provider_metadata_are_used(self, ubiquitin_pdb, tmp_settings):
        # The coordinate file cross-references P62988; the RCSB polymer entity
        # additionally reports P0CG48. Either should satisfy the query.
        report = validate_structure(
            ubiquitin_pdb,
            query_sequence=UBIQUITIN,
            query_accession="P0CG48",
            settings=tmp_settings,
            extra_accessions=["P0CG48"],
        )
        assert report.mapping_level is MappingLevel.M4_EXACT

    def test_fragment_of_the_query_keeps_identity_but_loses_coverage(
        self, ubiquitin_pdb, tmp_settings
    ):
        longer_query = UBIQUITIN + "GGGGSGGGGS" * 8
        report = validate_structure(
            ubiquitin_pdb,
            query_sequence=longer_query,
            query_accession="P62988",
            settings=tmp_settings,
        )
        assert report.identity_percent == pytest.approx(100.0)
        assert report.coverage_percent < 100.0
        assert report.mapping_level is not MappingLevel.M4_EXACT

    def test_chain_and_lengths_are_reported(self, ubiquitin_pdb, tmp_settings):
        report = validate_structure(
            ubiquitin_pdb,
            query_sequence=UBIQUITIN,
            query_accession="P62988",
            settings=tmp_settings,
        )
        assert report.subject_chain == "A"
        assert report.query_length == 76
        assert report.subject_length == 76


class TestMessages:
    def test_no_result_report(self):
        report = no_result_report()
        assert report.mapping_level is MappingLevel.M0_NO_RESULT
        assert report.message

    def test_describe_exact_mentions_the_accession(self):
        message = describe(MappingLevel.M4_EXACT, 100.0, 100.0, True, "P24941", ["P24941"])
        assert "P24941" in message and "Exact" in message

    def test_describe_related_carries_a_warning(self):
        message = describe(MappingLevel.M2_SEQUENCE_MATCH, 62.0, 80.0, False, "P01308", ["P14735"])
        assert message.startswith("WARNING")
        assert "P14735" in message
