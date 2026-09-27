"""Structure providers and the fallback manager, driven by a fake transport."""

from __future__ import annotations

import json

import pytest

from bioseqinsight.models.results import (
    AttemptRecord,
    MappingLevel,
    RetrievalStatus,
    StructureResult,
)
from bioseqinsight.services.cache import ResponseCache
from bioseqinsight.services.transport import FakeTransport, ScriptedResponse
from bioseqinsight.structures.alphafold import AlphaFoldProvider
from bioseqinsight.structures.base import (
    classify_query,
    is_pdb_id,
    is_uniprot_accession,
    normalise_accession,
)
from bioseqinsight.structures.esmatlas import ESMAtlasProvider
from bioseqinsight.structures.manager import StructureManager
from bioseqinsight.structures.rcsb import RCSBProvider
from bioseqinsight.structures.uniprot import UniProtClient
from conftest import UBIQUITIN, make_client

SEARCH_HIT = json.dumps({"result_set": [{"identifier": "1UBQ_1", "score": 1.0}]})
ENTITY = json.dumps(
    {
        "rcsb_polymer_entity_container_identifiers": {
            "reference_sequence_identifiers": [
                {"database_name": "UniProt", "database_accession": "P0CG48"}
            ]
        }
    }
)


class TestQueryClassification:
    @pytest.mark.parametrize("value", ["P24941", "P69905", "Q30597", "A0A0B4J2F0"])
    def test_uniprot_accessions_are_recognised(self, value):
        assert is_uniprot_accession(value)
        assert classify_query(value) == "uniprot"

    @pytest.mark.parametrize("value", ["1UBQ", "1aq1", "4HHB"])
    def test_pdb_identifiers_are_recognised(self, value):
        assert is_pdb_id(value)
        assert classify_query(value) == "pdb_id"

    def test_a_sequence_is_recognised(self):
        assert classify_query(UBIQUITIN) == "sequence"

    def test_nonsense_is_unknown(self):
        assert classify_query("???") == "unknown"
        assert classify_query("") == "unknown"

    @pytest.mark.parametrize(
        "raw,expected",
        [("sp|P24941|CDK2_HUMAN", "P24941"), ("AF-P69905-F1", "P69905"), (" p24941 ", "P24941")],
    )
    def test_accession_normalisation(self, raw, expected):
        assert normalise_accession(raw) == expected


class TestAlphaFoldProvider:
    """The provider queries the prediction API first, then the file URL it
    returns. This is what fixed the real defect where hard-coded model
    versions (v4, v3) all 404'd once EBI moved the database past them: every
    accession looked like a coverage gap instead of the version guess being
    wrong. See the module docstring in structures/alphafold.py.
    """

    def _provider(self, transport, settings):
        return AlphaFoldProvider(make_client(transport), settings, None)

    def _api_payload(self, accession: str, pdb_url: str, entry_id: str | None = None) -> str:
        return json.dumps(
            [{"entryId": entry_id or f"AF-{accession}-F1", "uniprotAccession": accession,
              "pdbUrl": pdb_url, "latestVersion": 4}]
        )

    def test_successful_retrieval(self, tmp_settings, ubiquitin_pdb):
        pdb_url = "https://alphafold.ebi.ac.uk/files/AF-P0CG48-F1-model_v4.pdb"
        transport = FakeTransport(
            {
                "api/prediction/P0CG48": [
                    ScriptedResponse(200, self._api_payload("P0CG48", pdb_url))
                ],
                pdb_url: [ScriptedResponse(200, ubiquitin_pdb)],
            }
        )
        result = self._provider(transport, tmp_settings).retrieve("P0CG48", accession="P0CG48")
        assert result.status is RetrievalStatus.OK
        assert result.confidence_kind == "plddt"
        assert result.download_path
        assert result.identifier == "AF-P0CG48-F1"

    def test_404_from_the_prediction_api_is_a_coverage_gap_not_a_service_failure(
        self, tmp_settings
    ):
        transport = FakeTransport(
            {"api/prediction/P00720": [ScriptedResponse(404, "no model")]}
        )
        result = self._provider(transport, tmp_settings).retrieve("P00720", accession="P00720")
        assert result.status is RetrievalStatus.NOT_FOUND
        assert "coverage gap" in result.error

    def test_empty_prediction_list_is_a_coverage_gap(self, tmp_settings):
        transport = FakeTransport(
            {"api/prediction/P00720": [ScriptedResponse(200, "[]")]}
        )
        result = self._provider(transport, tmp_settings).retrieve("P00720", accession="P00720")
        assert result.status is RetrievalStatus.NOT_FOUND

    def test_api_does_not_hard_code_a_model_version(self, tmp_settings, ubiquitin_pdb):
        # The whole point of the API-first design: whatever version the API
        # names in pdbUrl is fetched, with no guessing on this side.
        pdb_url = "https://alphafold.ebi.ac.uk/files/AF-P0CG48-F1-model_v9.pdb"
        transport = FakeTransport(
            {
                "api/prediction/P0CG48": [
                    ScriptedResponse(200, self._api_payload("P0CG48", pdb_url))
                ],
                pdb_url: [ScriptedResponse(200, ubiquitin_pdb)],
            }
        )
        result = self._provider(transport, tmp_settings).retrieve("P0CG48", accession="P0CG48")
        assert result.status is RetrievalStatus.OK
        assert result.structure_url == pdb_url

    def test_a_clean_two_call_success_is_not_reported_as_recovered(
        self, tmp_settings, ubiquitin_pdb
    ):
        # Regression test for a real defect: the API-discovery-then-download
        # design means a completely clean run makes two HTTP calls, and
        # StructureResult.recovered must not confuse "two calls" with "a
        # retry after failure". A live benchmark run once showed 67/67
        # successful AlphaFold retrievals mislabelled as recovered because
        # of exactly this.
        pdb_url = "https://alphafold.ebi.ac.uk/files/AF-P0CG48-F1-model_v4.pdb"
        transport = FakeTransport(
            {
                "api/prediction/P0CG48": [
                    ScriptedResponse(200, self._api_payload("P0CG48", pdb_url))
                ],
                pdb_url: [ScriptedResponse(200, ubiquitin_pdb)],
            }
        )
        result = self._provider(transport, tmp_settings).retrieve("P0CG48", accession="P0CG48")
        assert result.status is RetrievalStatus.OK
        assert len(result.attempts) == 2  # both calls were needed and both succeeded
        assert result.recovered is False

    def test_a_failed_api_call_that_then_succeeds_is_reported_as_recovered(
        self, tmp_settings, ubiquitin_pdb
    ):
        pdb_url = "https://alphafold.ebi.ac.uk/files/AF-P0CG48-F1-model_v4.pdb"
        transport = FakeTransport(
            {
                "api/prediction/P0CG48": [
                    ScriptedResponse(503, "down"),
                    ScriptedResponse(200, self._api_payload("P0CG48", pdb_url)),
                ],
                pdb_url: [ScriptedResponse(200, ubiquitin_pdb)],
            }
        )
        result = self._provider(transport, tmp_settings).retrieve("P0CG48", accession="P0CG48")
        assert result.status is RetrievalStatus.OK
        assert result.recovered is True

    def test_malformed_api_response_is_a_service_failure(self, tmp_settings):
        transport = FakeTransport(
            {"api/prediction/P0CG48": [ScriptedResponse(200, "not json")]}
        )
        result = self._provider(transport, tmp_settings).retrieve("P0CG48", accession="P0CG48")
        assert result.status is RetrievalStatus.SERVICE_ERROR

    def test_file_download_failing_after_a_successful_api_lookup_is_a_service_failure(
        self, tmp_settings
    ):
        pdb_url = "https://alphafold.ebi.ac.uk/files/AF-P0CG48-F1-model_v4.pdb"
        transport = FakeTransport(
            {
                "api/prediction/P0CG48": [
                    ScriptedResponse(200, self._api_payload("P0CG48", pdb_url))
                ],
                pdb_url: [ScriptedResponse(503, "down")],
            }
        )
        result = self._provider(transport, tmp_settings).retrieve("P0CG48", accession="P0CG48")
        assert result.status is RetrievalStatus.SERVICE_ERROR

    def test_missing_accession_is_a_user_error(self, tmp_settings):
        result = self._provider(FakeTransport({}), tmp_settings).retrieve("", accession="")
        assert result.status is RetrievalStatus.INVALID_INPUT

    def test_offline_mode_skips_the_request(self, tmp_settings):
        tmp_settings.offline = True
        transport = FakeTransport({})
        result = self._provider(transport, tmp_settings).retrieve("P0CG48", accession="P0CG48")
        assert result.status is RetrievalStatus.SKIPPED
        assert transport.calls == []

    def test_supports_only_accessions(self, tmp_settings):
        provider = self._provider(FakeTransport({}), tmp_settings)
        assert provider.supports("uniprot") and not provider.supports("sequence")


class TestESMAtlasProvider:
    def _provider(self, transport, settings):
        return ESMAtlasProvider(make_client(transport), settings, None)

    def test_successful_fold(self, tmp_settings, ubiquitin_pdb):
        transport = FakeTransport({"esmatlas.com": [ScriptedResponse(200, ubiquitin_pdb)]})
        result = self._provider(transport, tmp_settings).retrieve(UBIQUITIN, sequence=UBIQUITIN)
        assert result.status is RetrievalStatus.OK
        assert result.confidence_kind == "plddt"

    def test_repeated_504_is_reported_as_a_service_error_with_attempts(self, tmp_settings):
        transport = FakeTransport({"esmatlas.com": [ScriptedResponse(504, "gateway timeout")]})
        result = self._provider(transport, tmp_settings).retrieve(UBIQUITIN, sequence=UBIQUITIN)
        assert result.status is RetrievalStatus.SERVICE_ERROR
        assert len(result.attempts) == tmp_settings.max_attempts
        assert "504" in result.error

    def test_recovery_after_transient_504(self, tmp_settings, ubiquitin_pdb):
        transport = FakeTransport(
            {
                "esmatlas.com": [
                    ScriptedResponse(504, "gateway timeout"),
                    ScriptedResponse(200, ubiquitin_pdb),
                ]
            }
        )
        result = self._provider(transport, tmp_settings).retrieve(UBIQUITIN, sequence=UBIQUITIN)
        assert result.status is RetrievalStatus.OK
        assert result.recovered is True

    def test_oversized_sequence_is_refused_before_the_request(self, tmp_settings):
        tmp_settings.esm_max_length = 50
        transport = FakeTransport({})
        result = self._provider(transport, tmp_settings).retrieve(UBIQUITIN, sequence=UBIQUITIN)
        assert result.status is RetrievalStatus.INVALID_INPUT
        assert transport.calls == []

    def test_invalid_sequence_is_rejected(self, tmp_settings):
        result = self._provider(FakeTransport({}), tmp_settings).retrieve("", sequence="")
        assert result.status is RetrievalStatus.INVALID_INPUT

    def test_html_error_page_is_not_accepted_as_a_structure(self, tmp_settings):
        transport = FakeTransport({"esmatlas.com": [ScriptedResponse(200, "<html>oops</html>")]})
        result = self._provider(transport, tmp_settings).retrieve(UBIQUITIN, sequence=UBIQUITIN)
        assert result.status is not RetrievalStatus.OK


class TestRCSBProvider:
    def _provider(self, transport, settings, cache=None):
        return RCSBProvider(make_client(transport), settings, cache)

    def test_direct_entry_download(self, tmp_settings, ubiquitin_pdb):
        transport = FakeTransport({"files.rcsb.org": [ScriptedResponse(200, ubiquitin_pdb)]})
        result = self._provider(transport, tmp_settings).fetch_entry("1UBQ")
        assert result.status is RetrievalStatus.OK
        assert result.identifier == "1UBQ"
        assert result.experimental is True
        assert result.confidence_kind == "bfactor"

    def test_malformed_identifier_is_rejected(self, tmp_settings):
        result = self._provider(FakeTransport({}), tmp_settings).fetch_entry("NOPE")
        assert result.status is RetrievalStatus.INVALID_INPUT

    def test_missing_entry_is_not_found(self, tmp_settings):
        transport = FakeTransport({"files.rcsb.org": [ScriptedResponse(404, "missing")]})
        assert self._provider(transport, tmp_settings).fetch_entry("9ZZZ").status is (
            RetrievalStatus.NOT_FOUND
        )

    def test_sequence_search_downloads_the_hit(self, tmp_settings, ubiquitin_pdb):
        transport = FakeTransport(
            {
                "search.rcsb.org": [ScriptedResponse(200, SEARCH_HIT)],
                "data.rcsb.org": [ScriptedResponse(200, ENTITY)],
                "files.rcsb.org": [ScriptedResponse(200, ubiquitin_pdb)],
            }
        )
        result = self._provider(transport, tmp_settings).search_by_sequence(
            UBIQUITIN, accession="P0CG48"
        )
        assert result.status is RetrievalStatus.OK
        assert result.identifier == "1UBQ"
        assert "P0CG48" in result.cross_references

    def test_empty_search_result_is_not_found(self, tmp_settings):
        transport = FakeTransport(
            {"search.rcsb.org": [ScriptedResponse(200, json.dumps({"result_set": []}))]}
        )
        result = self._provider(transport, tmp_settings).search_by_sequence(UBIQUITIN)
        assert result.status is RetrievalStatus.NOT_FOUND

    def test_short_sequence_is_refused(self, tmp_settings):
        result = self._provider(FakeTransport({}), tmp_settings).search_by_sequence("MKT")
        assert result.status is RetrievalStatus.INVALID_INPUT

    def test_search_prefers_a_hit_carrying_the_queried_accession(self, tmp_settings, ubiquitin_pdb):
        two_hits = json.dumps(
            {
                "result_set": [
                    {"identifier": "9XXX_1", "score": 1.0},
                    {"identifier": "1UBQ_1", "score": 0.5},
                ]
            }
        )
        transport = FakeTransport(
            {
                "search.rcsb.org": [ScriptedResponse(200, two_hits)],
                "polymer_entity/9XXX": [
                    ScriptedResponse(
                        200,
                        json.dumps(
                            {
                                "rcsb_polymer_entity_container_identifiers": {
                                    "reference_sequence_identifiers": [
                                        {
                                            "database_name": "UniProt",
                                            "database_accession": "P99999",
                                        }
                                    ]
                                }
                            }
                        ),
                    )
                ],
                "polymer_entity/1UBQ": [ScriptedResponse(200, ENTITY)],
                "files.rcsb.org": [ScriptedResponse(200, ubiquitin_pdb)],
            }
        )
        result = self._provider(transport, tmp_settings).search_by_sequence(
            UBIQUITIN, accession="P0CG48"
        )
        # The top-scoring hit is 9XXX, but it is not the requested accession.
        assert result.identifier == "1UBQ"

    def test_entity_lookup_failure_is_not_fatal(self, tmp_settings, ubiquitin_pdb):
        transport = FakeTransport(
            {
                "search.rcsb.org": [ScriptedResponse(200, SEARCH_HIT)],
                "data.rcsb.org": [ScriptedResponse(500, "boom")],
                "files.rcsb.org": [ScriptedResponse(200, ubiquitin_pdb)],
            }
        )
        result = self._provider(transport, tmp_settings).search_by_sequence(UBIQUITIN)
        assert result.status is RetrievalStatus.OK
        assert result.cross_references == []


class TestUniProtClient:
    def test_sequence_is_fetched_and_cleaned(self, tmp_settings):
        transport = FakeTransport(
            {"rest.uniprot.org": [ScriptedResponse(200, f">sp|P0CG48|UBC\n{UBIQUITIN}")]}
        )
        client = UniProtClient(make_client(transport), tmp_settings, None)
        assert client.fetch_sequence("P0CG48") == UBIQUITIN

    def test_missing_entry_raises_not_found(self, tmp_settings):
        from bioseqinsight.services.errors import NotFoundError

        transport = FakeTransport({"rest.uniprot.org": [ScriptedResponse(404, "nope")]})
        client = UniProtClient(make_client(transport), tmp_settings, None)
        with pytest.raises(NotFoundError):
            client.fetch_sequence("P00000")

    def test_entry_metadata_is_parsed(self, tmp_settings):
        payload = json.dumps(
            {
                "primaryAccession": "P0CG48",
                "sequence": {"value": UBIQUITIN},
                "organism": {"scientificName": "Homo sapiens"},
                "proteinDescription": {"recommendedName": {"fullName": {"value": "Polyubiquitin-C"}}},
                "uniProtKBCrossReferences": [{"database": "PDB", "id": "1UBQ"}],
            }
        )
        transport = FakeTransport({"rest.uniprot.org": [ScriptedResponse(200, payload)]})
        client = UniProtClient(make_client(transport), tmp_settings, None)
        entry = client.fetch_entry("P0CG48")
        assert entry.organism == "Homo sapiens"
        assert entry.name == "Polyubiquitin-C"
        assert entry.pdb_ids == ["1UBQ"]
        assert entry.length == 76


class TestStructureManager:
    def _manager(self, transport, settings):
        return StructureManager(
            settings, make_client(transport), ResponseCache(settings.cache_dir)
        )

    def test_accession_query_resolves_to_an_exact_match(self, tmp_settings, working_transport):
        outcome = self._manager(working_transport, tmp_settings).retrieve("P0CG48")
        assert outcome.best_mapping is MappingLevel.M4_EXACT
        assert outcome.warnings == []

    def test_exact_match_stops_the_fallback_chain(self, tmp_settings, working_transport):
        outcome = self._manager(working_transport, tmp_settings).retrieve("P0CG48")
        assert len(outcome.results) == 1  # rcsb is first and answered exactly

    def test_stop_on_exact_can_be_disabled(self, tmp_settings, working_transport):
        outcome = self._manager(working_transport, tmp_settings).retrieve(
            "P0CG48", stop_on_exact=False
        )
        assert {r.source for r in outcome.results} == {"rcsb", "alphafold", "esmatlas"}

    def test_fallback_when_the_first_provider_fails(self, tmp_settings, ubiquitin_pdb):
        pdb_url = "https://alphafold.ebi.ac.uk/files/AF-P0CG48-F1-model_v4.pdb"
        transport = FakeTransport(
            {
                "rest.uniprot.org": [ScriptedResponse(200, f">sp|P0CG48|UBC\n{UBIQUITIN}")],
                "search.rcsb.org": [ScriptedResponse(503, "down")],
                "api/prediction": [
                    ScriptedResponse(
                        200,
                        json.dumps(
                            [{"entryId": "AF-P0CG48-F1", "pdbUrl": pdb_url, "latestVersion": 4}]
                        ),
                    )
                ],
                pdb_url: [ScriptedResponse(200, ubiquitin_pdb)],
            }
        )
        outcome = self._manager(transport, tmp_settings).retrieve("P0CG48")
        sources = [r.source for r in outcome.results]
        assert sources[0] == "rcsb" and not outcome.results[0].ok
        assert "alphafold" in sources
        assert outcome.best.source == "alphafold"

    def test_all_providers_failing_gives_m0_and_a_warning(self, tmp_settings):
        transport = FakeTransport(
            {
                "rest.uniprot.org": [ScriptedResponse(200, f">sp|P0CG48|UBC\n{UBIQUITIN}")],
                "search.rcsb.org": [ScriptedResponse(503, "down")],
                "alphafold.ebi.ac.uk": [ScriptedResponse(503, "down")],
                "esmatlas.com": [ScriptedResponse(504, "gateway timeout")],
            }
        )
        outcome = self._manager(transport, tmp_settings).retrieve("P0CG48")
        assert outcome.best is None
        assert outcome.best_mapping is MappingLevel.M0_NO_RESULT
        assert any("No resource returned" in w for w in outcome.warnings)

    def test_unrecognised_input_is_reported_without_any_request(self, tmp_settings):
        transport = FakeTransport({})
        outcome = self._manager(transport, tmp_settings).retrieve("???")
        assert outcome.results[0].status is RetrievalStatus.INVALID_INPUT
        assert transport.calls == []

    def test_pdb_identifier_query(self, tmp_settings, ubiquitin_pdb):
        transport = FakeTransport({"files.rcsb.org": [ScriptedResponse(200, ubiquitin_pdb)]})
        outcome = self._manager(transport, tmp_settings).retrieve("1UBQ")
        assert outcome.results[0].identifier == "1UBQ"
        # No query sequence is available for a bare PDB id, so identity is
        # unverified rather than silently assumed.
        assert outcome.best_mapping is MappingLevel.M1_UNVERIFIED

    def test_failed_uniprot_lookup_degrades_gracefully(self, tmp_settings, ubiquitin_pdb):
        transport = FakeTransport(
            {
                "rest.uniprot.org": [ScriptedResponse(500, "down")],
                "alphafold.ebi.ac.uk": [ScriptedResponse(200, ubiquitin_pdb)],
            }
        )
        outcome = self._manager(transport, tmp_settings).retrieve("P0CG48")
        assert any("canonical sequence" in note for note in outcome.warnings)

    def test_second_identical_query_is_served_from_cache(self, tmp_settings, working_transport):
        manager = self._manager(working_transport, tmp_settings)
        manager.retrieve("P0CG48")
        calls_after_first = len(working_transport.calls)
        outcome = manager.retrieve("P0CG48")
        assert len(working_transport.calls) == calls_after_first
        assert outcome.best.from_cache is True

    def test_refresh_bypasses_the_cache(self, tmp_settings, working_transport):
        manager = self._manager(working_transport, tmp_settings)
        manager.retrieve("P0CG48")
        calls_after_first = len(working_transport.calls)
        manager.retrieve("P0CG48", refresh=True)
        assert len(working_transport.calls) > calls_after_first

    def test_offline_mode_makes_no_requests(self, tmp_settings):
        tmp_settings.offline = True
        transport = FakeTransport({})
        outcome = self._manager(transport, tmp_settings).retrieve(UBIQUITIN)
        assert transport.calls == []
        assert all(r.status is RetrievalStatus.SKIPPED for r in outcome.results)

    def test_outcome_is_json_serialisable(self, tmp_settings, working_transport):
        outcome = self._manager(working_transport, tmp_settings).retrieve("P0CG48")
        assert json.loads(json.dumps(outcome.to_dict()))["best_mapping"] == "M4"


class TestStructureResultRecovered:
    """StructureResult.recovered in isolation from any provider.

    A regression suite for a real defect: a live benchmark run once reported
    67 of 67 successful AlphaFold retrievals as "recovered after failure"
    when none of them had actually failed. The cause was conflating "more
    than one HTTP call was made" with "a call failed and was retried" --
    true for single-call providers by construction, false for AlphaFold DB,
    which always makes two calls (discover the model URL, then download it)
    even on a perfectly clean run.
    """

    def _attempt(self, attempt=1, status_code=200, error=None) -> AttemptRecord:
        return AttemptRecord(
            attempt=attempt, url="https://example.org", status_code=status_code,
            error=error, elapsed_s=0.01,
        )

    def _result(self, attempts, status=RetrievalStatus.OK) -> StructureResult:
        return StructureResult(source="test", query="q", query_type="sequence",
                                status=status, attempts=attempts)

    def test_a_single_successful_attempt_is_not_recovered(self):
        assert self._result([self._attempt()]).recovered is False

    def test_two_successful_attempts_are_not_recovered(self):
        # This is the exact case that was wrong: two clean calls, e.g. an
        # API discovery call followed by a file download, neither of which
        # failed.
        assert self._result([self._attempt(1), self._attempt(1)]).recovered is False

    def test_a_failed_attempt_followed_by_success_is_recovered(self):
        attempts = [self._attempt(1, status_code=503, error="service unavailable"),
                    self._attempt(2, status_code=200)]
        assert self._result(attempts).recovered is True

    def test_an_error_without_a_status_code_still_counts_as_a_failure(self):
        attempts = [self._attempt(1, status_code=None, error="timeout after 5s"),
                    self._attempt(2, status_code=200)]
        assert self._result(attempts).recovered is True

    def test_a_failed_result_is_never_recovered_even_with_failed_attempts(self):
        attempts = [self._attempt(1, status_code=503, error="down")]
        assert self._result(attempts, status=RetrievalStatus.SERVICE_ERROR).recovered is False

    def test_three_clean_calls_are_still_not_recovered(self):
        # Guards against a future provider that legitimately makes more than
        # two calls on a clean path.
        attempts = [self._attempt(1), self._attempt(1), self._attempt(1)]
        assert self._result(attempts).recovered is False
