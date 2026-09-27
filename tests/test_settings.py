"""Configuration loading and validation."""

from __future__ import annotations

import json
import os

import pytest

from bioseqinsight.config.settings import Settings


class TestDefaults:
    def test_defaults_are_valid(self):
        Settings().validate()

    def test_thresholds_are_ordered_by_default(self):
        settings = Settings()
        assert (
            settings.related_identity_threshold
            <= settings.high_identity_threshold
            <= settings.exact_identity_threshold
        )

    def test_cache_ttl_converts_days_to_seconds(self):
        assert Settings(cache_ttl_days=2).cache_ttl_s == 2 * 24 * 3600


class TestValidation:
    def test_zero_attempts_is_rejected(self):
        with pytest.raises(ValueError, match="max_attempts"):
            Settings(max_attempts=0).validate()

    def test_negative_timeout_is_rejected(self):
        with pytest.raises(ValueError, match="timeouts"):
            Settings(http_timeout_s=-1).validate()

    def test_identity_outside_a_percentage_is_rejected(self):
        with pytest.raises(ValueError, match="percentage"):
            Settings(exact_identity_threshold=150).validate()

    def test_inverted_thresholds_are_rejected(self):
        with pytest.raises(ValueError, match="high_identity_threshold"):
            Settings(high_identity_threshold=100, exact_identity_threshold=90).validate()

    def test_unknown_provider_is_rejected(self):
        with pytest.raises(ValueError, match="Unknown providers"):
            Settings(fallback_order=("rcsb", "magic")).validate()

    def test_zero_workers_is_rejected(self):
        with pytest.raises(ValueError, match="batch_workers"):
            Settings(batch_workers=0).validate()


class TestLoading:
    def test_from_dict_ignores_unknown_keys(self):
        settings = Settings.from_dict({"max_attempts": 5, "nonsense": True})
        assert settings.max_attempts == 5

    def test_string_values_are_coerced(self):
        settings = Settings.from_dict(
            {"max_attempts": "7", "http_timeout_s": "12.5", "offline": "true"}
        )
        assert settings.max_attempts == 7
        assert settings.http_timeout_s == 12.5
        assert settings.offline is True

    def test_comma_separated_tuples_are_parsed(self):
        settings = Settings.from_dict({"fallback_order": "alphafold,rcsb"})
        assert settings.fallback_order == ("alphafold", "rcsb")

    def test_file_is_loaded(self, tmp_path):
        path = tmp_path / "settings.json"
        path.write_text(json.dumps({"max_attempts": 9, "offline": True}))
        settings = Settings.load(path=path, use_env=False)
        assert settings.max_attempts == 9
        assert settings.source == str(path)

    def test_malformed_file_is_reported(self, tmp_path):
        path = tmp_path / "settings.json"
        path.write_text("{not json")
        with pytest.raises(ValueError, match="Cannot read settings file"):
            Settings.load(path=path, use_env=False)

    def test_environment_overrides_the_file(self, tmp_path):
        path = tmp_path / "settings.json"
        path.write_text(json.dumps({"max_attempts": 2}))
        os.environ["BIOSEQINSIGHT_MAX_ATTEMPTS"] = "8"
        try:
            assert Settings.load(path=path).max_attempts == 8
        finally:
            del os.environ["BIOSEQINSIGHT_MAX_ATTEMPTS"]

    def test_explicit_overrides_win(self, tmp_path):
        path = tmp_path / "settings.json"
        path.write_text(json.dumps({"max_attempts": 2}))
        assert Settings.load(path=path, use_env=False, max_attempts=4).max_attempts == 4

    def test_none_overrides_are_ignored(self):
        assert Settings.load(use_env=False, max_attempts=None).max_attempts == 3


class TestPersistence:
    def test_save_then_load(self, tmp_path):
        original = Settings(max_attempts=6, offline=True)
        path = original.save(tmp_path / "s.json")
        reloaded = Settings.load(path=path, use_env=False)
        assert reloaded.max_attempts == 6 and reloaded.offline is True

    def test_to_dict_is_json_serialisable(self):
        assert json.loads(json.dumps(Settings().to_dict()))["max_attempts"] == 3

    def test_tuples_become_lists_for_json(self):
        assert isinstance(Settings().to_dict()["fallback_order"], list)

    def test_ensure_directories_creates_them(self, tmp_path):
        settings = Settings(
            cache_dir=str(tmp_path / "c"),
            download_dir=str(tmp_path / "d"),
            projects_dir=str(tmp_path / "p"),
            log_file=str(tmp_path / "logs" / "a.log"),
        )
        settings.ensure_directories()
        assert (tmp_path / "c").is_dir() and (tmp_path / "logs").is_dir()
