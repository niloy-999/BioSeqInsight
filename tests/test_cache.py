"""On-disk response cache."""

from __future__ import annotations

import json
import time

import pytest

from bioseqinsight.services.cache import CacheEntry, NullCache, ResponseCache, make_key
from bioseqinsight.services.errors import CacheError


class TestKeys:
    def test_same_inputs_give_the_same_key(self):
        assert make_key("rcsb", "P24941", v=1) == make_key("rcsb", "P24941", v=1)

    def test_different_query_gives_a_different_key(self):
        assert make_key("rcsb", "P24941") != make_key("rcsb", "P69905")

    def test_different_namespace_gives_a_different_key(self):
        assert make_key("rcsb", "X") != make_key("alphafold", "X")

    def test_parameter_order_does_not_matter(self):
        assert make_key("n", "q", a=1, b=2) == make_key("n", "q", b=2, a=1)

    def test_parameters_affect_the_key(self):
        assert make_key("n", "q", cutoff=0.9) != make_key("n", "q", cutoff=0.5)


class TestRoundTrip:
    def test_store_then_retrieve(self, cache):
        key = make_key("rcsb", "1UBQ")
        cache.put("rcsb", key, "ATOM payload", query="1UBQ", source="rcsb", url="https://x")
        payload, entry = cache.get("rcsb", key)
        assert payload == "ATOM payload"
        assert entry.query == "1UBQ"
        assert entry.url == "https://x"
        assert entry.software_version

    def test_miss_returns_none(self, cache):
        assert cache.get("rcsb", "nope") is None

    def test_hit_and_miss_counters(self, cache):
        key = make_key("rcsb", "1UBQ")
        cache.put("rcsb", key, "x", query="1UBQ", source="rcsb")
        cache.get("rcsb", key)
        cache.get("rcsb", "absent")
        assert cache.hits == 1 and cache.misses == 1
        assert cache.stats()["hit_rate_percent"] == pytest.approx(50.0)

    def test_refresh_bypasses_a_live_entry(self, cache):
        key = make_key("rcsb", "1UBQ")
        cache.put("rcsb", key, "x", query="1UBQ", source="rcsb")
        assert cache.get("rcsb", key, refresh=True) is None

    def test_metadata_file_is_human_readable(self, cache):
        key = make_key("rcsb", "1UBQ")
        cache.put("rcsb", key, "payload", query="1UBQ", source="rcsb")
        meta = json.loads((cache.directory / "rcsb" / f"{key}.json").read_text())
        assert meta["query"] == "1UBQ"
        assert meta["size_bytes"] == len("payload")


class TestIntegrity:
    def test_tampered_payload_is_rejected(self, cache):
        key = make_key("rcsb", "1UBQ")
        cache.put("rcsb", key, "original", query="1UBQ", source="rcsb")
        (cache.directory / "rcsb" / f"{key}.payload").write_text("tampered")
        assert cache.get("rcsb", key) is None

    def test_corrupt_metadata_is_discarded(self, cache):
        key = make_key("rcsb", "1UBQ")
        cache.put("rcsb", key, "payload", query="1UBQ", source="rcsb")
        (cache.directory / "rcsb" / f"{key}.json").write_text("{not json")
        assert cache.get("rcsb", key) is None

    def test_expired_entries_are_not_served(self, tmp_path):
        cache = ResponseCache(tmp_path / "c", ttl_s=0.0)
        key = make_key("rcsb", "1UBQ")
        cache.put("rcsb", key, "payload", query="1UBQ", source="rcsb")
        time.sleep(0.01)
        assert cache.get("rcsb", key) is None

    def test_no_ttl_never_expires(self, tmp_path):
        cache = ResponseCache(tmp_path / "c", ttl_s=None)
        key = make_key("rcsb", "1UBQ")
        cache.put("rcsb", key, "payload", query="1UBQ", source="rcsb")
        assert cache.get("rcsb", key) is not None


class TestManagement:
    def test_entries_are_listed(self, cache):
        for accession in ("P1", "P2"):
            cache.put("rcsb", make_key("rcsb", accession), "x", query=accession, source="rcsb")
        assert len(cache.entries("rcsb")) == 2
        assert len(cache.entries()) == 2

    def test_delete_removes_both_files(self, cache):
        key = make_key("rcsb", "1UBQ")
        cache.put("rcsb", key, "x", query="1UBQ", source="rcsb")
        assert cache.delete("rcsb", key) is True
        assert cache.get("rcsb", key) is None

    def test_deleting_a_missing_entry_is_harmless(self, cache):
        assert cache.delete("rcsb", "absent") is False

    def test_clear_a_single_namespace(self, cache):
        cache.put("rcsb", make_key("rcsb", "a"), "x", query="a", source="rcsb")
        cache.put("alphafold", make_key("alphafold", "b"), "y", query="b", source="af")
        assert cache.clear("rcsb") == 1
        assert len(cache.entries()) == 1

    def test_clear_everything(self, cache):
        cache.put("rcsb", make_key("rcsb", "a"), "x", query="a", source="rcsb")
        cache.put("alphafold", make_key("alphafold", "b"), "y", query="b", source="af")
        assert cache.clear() == 2
        assert cache.entries() == []


class TestDisabled:
    def test_disabled_cache_never_returns_a_hit(self, tmp_path):
        cache = ResponseCache(tmp_path / "c", enabled=False)
        with pytest.raises(CacheError, match="disabled"):
            cache.put("rcsb", "k", "x", query="q", source="rcsb")
        assert cache.get("rcsb", "k") is None

    def test_null_cache_is_inert(self):
        null = NullCache()
        assert null.get("any", "key") is None
        entry = null.put("any", key="key", payload="x", query="q", source="s")
        assert isinstance(entry, CacheEntry)
