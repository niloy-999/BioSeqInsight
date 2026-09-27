"""GUI output formatting.

The view classes need Tkinter, but the functions that turn a result object
into displayed text do not. They are pure functions for exactly this reason:
the wording the user sees, including the warnings that stop a non-exact match
being read as a confirmed one, is covered by tests that run headless in CI.
"""

from __future__ import annotations

import json

from bioseqinsight.core.protein import ProteinAnalyzer
from bioseqinsight.core.sequence import SequenceAnalyzer
from bioseqinsight.core.translation import six_frame_translation
from bioseqinsight.gui.formatting import (
    format_hydropathy,
    format_protein_result,
    format_sequence_result,
    format_six_frames,
    format_structure_outcome,
)
from bioseqinsight.services.cache import ResponseCache
from bioseqinsight.services.transport import FakeTransport, ScriptedResponse
from bioseqinsight.structures.manager import StructureManager
from conftest import UBIQUITIN, make_client

DNA = "ATG" + "GCTAGCTTACCGATGGC" * 4 + "TAA"


class TestSequenceFormatting:
    def test_key_fields_are_present(self):
        text = format_sequence_result(SequenceAnalyzer(DNA, "demo").analyze(min_orf_aa=5))
        assert "demo" in text
        assert "GC content" in text
        assert "Melting temp" in text
        assert "ORFs" in text

    def test_warnings_are_shown(self):
        text = format_sequence_result(SequenceAnalyzer("ATGCNNNNATGC").analyze())
        assert "Warnings" in text and "ambiguity" in text

    def test_tm_method_and_conditions_are_stated(self):
        text = format_sequence_result(SequenceAnalyzer(DNA).analyze())
        assert "nearest_neighbor" in text
        assert "mM Na+" in text

    def test_motif_hits_are_listed(self):
        result = SequenceAnalyzer(DNA).analyze(motif="GCTAGC")
        assert "Motif hits" in format_sequence_result(result)

    def test_output_ends_with_provenance(self):
        assert "computed locally" in format_sequence_result(SequenceAnalyzer(DNA).analyze())


class TestProteinFormatting:
    def test_key_fields_are_present(self):
        text = format_protein_result(ProteinAnalyzer(UBIQUITIN, "ubq").analyze())
        assert "Molecular weight" in text
        assert "Isoelectric point" in text
        assert "Extinction 280 nm" in text

    def test_sketch_carries_its_disclaimer(self):
        text = format_protein_result(ProteinAnalyzer(UBIQUITIN).analyze(include_sketch=True))
        assert "NOT a structure prediction" in text

    def test_composition_is_rendered(self):
        assert "Amino-acid composition" in format_protein_result(
            ProteinAnalyzer(UBIQUITIN).analyze()
        )

    def test_hydropathy_profile_renders_a_bar_chart(self):
        from bioseqinsight.core.protein import hydropathy_profile

        text = format_hydropathy(hydropathy_profile(UBIQUITIN, 9), 9)
        assert "hydropathy profile" in text
        assert "#" in text

    def test_short_sequence_hydropathy_message(self):
        assert "shorter than" in format_hydropathy([], 9)

    def test_six_frame_output_lists_all_frames(self):
        text = format_six_frames(six_frame_translation(DNA))
        for frame in ("+1", "+2", "+3", "-1", "-2", "-3"):
            assert f"Frame {frame}" in text


class TestStructureFormatting:
    def _outcome(self, transport, settings):
        manager = StructureManager(
            settings, make_client(transport), ResponseCache(settings.cache_dir)
        )
        return manager.retrieve("P0CG48")

    def test_exact_match_is_stated_plainly(self, tmp_settings, working_transport):
        text = format_structure_outcome(self._outcome(working_transport, tmp_settings))
        assert "MATCH STATUS" in text
        assert "M4" in text
        assert "Exact identity    YES" in text

    def test_attempts_and_recovery_are_shown(self, tmp_settings, ubiquitin_pdb):
        transport = FakeTransport(
            {
                "rest.uniprot.org": [ScriptedResponse(200, f">sp|P0CG48|UBC\n{UBIQUITIN}")],
                "search.rcsb.org": [
                    ScriptedResponse(503, "down"),
                    ScriptedResponse(
                        200, json.dumps({"result_set": [{"identifier": "1UBQ_1", "score": 1.0}]})
                    ),
                ],
                "data.rcsb.org": [ScriptedResponse(200, "{}")],
                "files.rcsb.org": [ScriptedResponse(200, ubiquitin_pdb)],
            }
        )
        text = format_structure_outcome(self._outcome(transport, tmp_settings))
        assert "attempts" in text

    def test_failure_is_reported_as_m0_with_warnings(self, tmp_settings):
        transport = FakeTransport(
            {
                "rest.uniprot.org": [ScriptedResponse(200, f">sp|P0CG48|UBC\n{UBIQUITIN}")],
                "search.rcsb.org": [ScriptedResponse(503, "down")],
                "alphafold.ebi.ac.uk": [ScriptedResponse(503, "down")],
                "esmatlas.com": [ScriptedResponse(504, "gateway timeout")],
            }
        )
        text = format_structure_outcome(self._outcome(transport, tmp_settings))
        assert "M0 - No result" in text
        assert "WARNINGS" in text

    def test_a_non_exact_match_is_never_presented_as_confirmed(self, tmp_settings, ubiquitin_pdb):
        # RCSB returns 1UBQ but reports a different accession: the interface
        # must say the identity was not confirmed.
        transport = FakeTransport(
            {
                "rest.uniprot.org": [ScriptedResponse(200, f">sp|P0CG48|UBC\n{UBIQUITIN}")],
                "search.rcsb.org": [
                    ScriptedResponse(
                        200, json.dumps({"result_set": [{"identifier": "1UBQ_1", "score": 1.0}]})
                    )
                ],
                "data.rcsb.org": [
                    ScriptedResponse(
                        200,
                        json.dumps(
                            {
                                "rcsb_polymer_entity_container_identifiers": {
                                    "reference_sequence_identifiers": [
                                        {"database_name": "UniProt", "database_accession": "P99999"}
                                    ]
                                }
                            }
                        ),
                    )
                ],
                "files.rcsb.org": [ScriptedResponse(200, ubiquitin_pdb)],
                "alphafold.ebi.ac.uk": [ScriptedResponse(404, "none")],
                "esmatlas.com": [ScriptedResponse(504, "timeout")],
            }
        )
        text = format_structure_outcome(self._outcome(transport, tmp_settings))
        assert "Exact identity    NO" in text
        assert "WARNINGS" in text
