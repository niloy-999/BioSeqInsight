"""Shared pytest fixtures.

Fixtures here guarantee that no test touches the network, the user's home
directory or the real cache. Every test in this suite is deterministic and
runs offline; that is a requirement, not an accident, because the CI matrix
must not depend on the availability of RCSB, EBI or ESM Atlas.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:  # allows `pytest` to run without installing first
    sys.path.insert(0, str(SRC))

from bioseqinsight.config.settings import Settings  # noqa: E402
from bioseqinsight.services.cache import ResponseCache  # noqa: E402
from bioseqinsight.services.http_client import HttpClient, RetryPolicy  # noqa: E402
from bioseqinsight.services.transport import FakeTransport, ScriptedResponse  # noqa: E402

DATA = Path(__file__).resolve().parent / "data"

UBIQUITIN = (
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
)
INSULIN_A = "GIVEQCCTSICSLYQLENYCN"
# A short synthetic construct with one unambiguous ORF on the forward strand.
DEMO_DNA = "TTTATGAAAGGGTTTCCCAAAGGGTTTCCCAAATAGCCC"


@pytest.fixture
def tmp_settings(tmp_path) -> Settings:
    """Settings confined to a temporary directory, with jitter disabled."""
    return Settings.from_dict(
        {
            "download_dir": str(tmp_path / "structures"),
            "cache_dir": str(tmp_path / "cache"),
            "projects_dir": str(tmp_path / "projects"),
            "log_file": str(tmp_path / "logs" / "test.log"),
            "cache_enabled": True,
            "jitter": False,
            "backoff_base_s": 0.0,
            "backoff_max_s": 0.0,
            "max_attempts": 3,
        },
        source="test-fixture",
    )


@pytest.fixture
def cache(tmp_settings) -> ResponseCache:
    return ResponseCache(tmp_settings.cache_dir, ttl_s=tmp_settings.cache_ttl_s)


def make_client(transport: FakeTransport, attempts: int = 3) -> HttpClient:
    """An HttpClient wired to a fake transport with instantaneous backoff."""
    return HttpClient(
        transport=transport,
        policy=RetryPolicy(
            max_attempts=attempts, backoff_base_s=0.0, backoff_max_s=0.0, jitter=False
        ),
        timeout_s=5.0,
        sleeper=lambda _seconds: None,
    )


@pytest.fixture
def ubiquitin_pdb() -> str:
    return (DATA / "1ubq_excerpt.pdb").read_text(encoding="utf-8")


@pytest.fixture
def working_transport(ubiquitin_pdb) -> FakeTransport:
    """A transport where every bundled service behaves correctly."""
    return FakeTransport(
        {
            "rest.uniprot.org": [ScriptedResponse(200, f">sp|P0CG48|UBC_HUMAN\n{UBIQUITIN}")],
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
                                    {"database_name": "UniProt", "database_accession": "P0CG48"}
                                ]
                            }
                        }
                    ),
                )
            ],
            "files.rcsb.org": [ScriptedResponse(200, ubiquitin_pdb)],
            # AlphaFold DB is queried in two steps: the prediction API names
            # the current model file, then that file is fetched. Both legs
            # need a scripted response.
            "api/prediction": [
                ScriptedResponse(
                    200,
                    json.dumps(
                        [
                            {
                                "entryId": "AF-P0CG48-F1",
                                "uniprotAccession": "P0CG48",
                                "pdbUrl": "https://alphafold.ebi.ac.uk/files/AF-P0CG48-F1-model_v4.pdb",
                                "latestVersion": 4,
                            }
                        ]
                    ),
                )
            ],
            "alphafold.ebi.ac.uk/files": [ScriptedResponse(200, ubiquitin_pdb)],
            "api.esmatlas.com": [ScriptedResponse(200, ubiquitin_pdb)],
        }
    )
