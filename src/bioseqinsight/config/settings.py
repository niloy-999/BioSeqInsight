"""Configuration management.

v1.0 hard-coded timeouts, output paths and endpoint URLs inside the modules
that used them. v2.0 resolves every one of those through :class:`Settings`,
which loads, in increasing order of precedence:

1. built-in defaults;
2. a JSON settings file (``~/.bioseqinsight/settings.json`` by default, or
   ``$BIOSEQINSIGHT_CONFIG``);
3. environment variables prefixed ``BIOSEQINSIGHT_``;
4. explicit keyword arguments (what the CLI flags and the GUI dialog set).

No credentials or API keys are involved: every bundled endpoint is a public,
anonymous HTTPS service.
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any

DEFAULT_CONFIG_DIR = Path.home() / ".bioseqinsight"
ENV_PREFIX = "BIOSEQINSIGHT_"

ESMATLAS_FOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
#: The prediction API is the authoritative source for the current model file
#: URL and is queried first; it returns whichever version is live without the
#: client having to know its number. ALPHAFOLD_FILE_URL is kept only as a
#: documented fallback template for anyone constructing a URL by hand; the
#: provider does not use it to discover models.
ALPHAFOLD_FILE_URL = "https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_{version}.pdb"
ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_POLYMER_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
UNIPROT_JSON_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"


@dataclass
class Settings:
    """Runtime configuration for one BioSeqInsight session."""

    # -- network ------------------------------------------------------------
    http_timeout_s: float = 60.0
    esm_timeout_s: float = 90.0
    max_attempts: int = 3
    backoff_base_s: float = 1.0
    backoff_max_s: float = 20.0
    jitter: bool = True
    offline: bool = False

    # -- cache --------------------------------------------------------------
    cache_enabled: bool = True
    cache_dir: str = str(DEFAULT_CONFIG_DIR / "cache")
    cache_ttl_days: float = 30.0

    # -- filesystem ---------------------------------------------------------
    download_dir: str = str(Path.cwd() / "structures")
    projects_dir: str = str(DEFAULT_CONFIG_DIR / "projects")
    log_file: str = str(DEFAULT_CONFIG_DIR / "logs" / "bioseqinsight.log")
    log_level: str = "INFO"
    log_json: bool = False

    # -- analysis defaults --------------------------------------------------
    genetic_code_table: int = 1
    min_orf_aa: int = 30
    tm_primer_nM: float = 250.0
    tm_na_mM: float = 50.0
    hydropathy_window: int = 9

    # -- structure retrieval ------------------------------------------------
    esm_max_length: int = 400
    #: Retained for backward compatibility with saved settings files and for
    #: anyone constructing an ALPHAFOLD_FILE_URL by hand. The AlphaFold
    #: provider itself discovers the current model version via
    #: ALPHAFOLD_API_URL and does not read this field.
    alphafold_versions: tuple[str, ...] = ("v4", "v3")
    fallback_order: tuple[str, ...] = ("rcsb", "alphafold", "esmatlas")
    rcsb_search_identity_cutoff: float = 0.9
    rcsb_search_max_hits: int = 5
    validate_structures: bool = True

    # -- mapping thresholds (the M0-M4 taxonomy) ---------------------------
    exact_identity_threshold: float = 99.0
    exact_coverage_threshold: float = 95.0
    high_identity_threshold: float = 95.0
    high_coverage_threshold: float = 80.0
    related_identity_threshold: float = 40.0
    related_coverage_threshold: float = 40.0

    # -- batch --------------------------------------------------------------
    batch_workers: int = 4
    batch_chunk_size: int = 100

    # -- visualisation ------------------------------------------------------
    viewer_command: str = ""  # empty means auto-detect
    open_browser: bool = True

    _source: str = field(default="defaults", repr=False)

    # -- construction -------------------------------------------------------

    @classmethod
    def load(
        cls,
        path: str | os.PathLike[str] | None = None,
        use_env: bool = True,
        **overrides: Any,
    ) -> Settings:
        """Build settings from file, environment and explicit overrides."""
        data: dict[str, Any] = {}
        source = "defaults"
        config_path = _resolve_config_path(path)
        if config_path and config_path.is_file():
            try:
                data.update(json.loads(config_path.read_text(encoding="utf-8")))
                source = str(config_path)
            except (OSError, json.JSONDecodeError) as exc:
                raise ValueError(f"Cannot read settings file {config_path}: {exc}") from exc
        if use_env:
            data.update(_from_environment())
            if any(k.startswith(ENV_PREFIX) for k in os.environ):
                source = f"{source}+env"
        data.update({k: v for k, v in overrides.items() if v is not None})
        return cls.from_dict(data, source=source)

    @classmethod
    def from_dict(cls, data: dict[str, Any], source: str = "dict") -> Settings:
        known = {f.name: f for f in fields(cls) if not f.name.startswith("_")}
        clean: dict[str, Any] = {}
        for key, value in data.items():
            if key not in known:
                continue
            clean[key] = _coerce(known[key].type, value)
        instance = cls(**clean)
        instance._source = source
        instance.validate()
        return instance

    # -- behaviour ----------------------------------------------------------

    def validate(self) -> None:
        if self.max_attempts < 1:
            raise ValueError("max_attempts must be at least 1")
        if self.http_timeout_s <= 0 or self.esm_timeout_s <= 0:
            raise ValueError("timeouts must be positive")
        if not 0 <= self.exact_identity_threshold <= 100:
            raise ValueError("exact_identity_threshold must be a percentage")
        if self.high_identity_threshold > self.exact_identity_threshold:
            raise ValueError("high_identity_threshold cannot exceed exact_identity_threshold")
        if self.related_identity_threshold > self.high_identity_threshold:
            raise ValueError("related_identity_threshold cannot exceed high_identity_threshold")
        if self.batch_workers < 1:
            raise ValueError("batch_workers must be at least 1")
        unknown = set(self.fallback_order) - {"rcsb", "alphafold", "esmatlas"}
        if unknown:
            raise ValueError(f"Unknown providers in fallback_order: {sorted(unknown)}")

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data.pop("_source", None)
        data["alphafold_versions"] = list(self.alphafold_versions)
        data["fallback_order"] = list(self.fallback_order)
        return data

    def save(self, path: str | os.PathLike[str] | None = None) -> Path:
        target = Path(path) if path else (DEFAULT_CONFIG_DIR / "settings.json")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True), encoding="utf-8")
        return target

    def ensure_directories(self) -> None:
        for directory in (self.cache_dir, self.download_dir, self.projects_dir):
            Path(directory).mkdir(parents=True, exist_ok=True)
        Path(self.log_file).parent.mkdir(parents=True, exist_ok=True)

    @property
    def cache_ttl_s(self) -> float:
        return self.cache_ttl_days * 24 * 3600

    @property
    def source(self) -> str:
        return self._source


def _resolve_config_path(path: str | os.PathLike[str] | None) -> Path | None:
    if path:
        return Path(path)
    env_path = os.environ.get(f"{ENV_PREFIX}CONFIG")
    if env_path:
        return Path(env_path)
    default = DEFAULT_CONFIG_DIR / "settings.json"
    return default if default.is_file() else None


def _from_environment() -> dict[str, Any]:
    out: dict[str, Any] = {}
    valid = {f.name for f in fields(Settings) if not f.name.startswith("_")}
    for key, value in os.environ.items():
        if not key.startswith(ENV_PREFIX):
            continue
        name = key[len(ENV_PREFIX) :].lower()
        if name in valid:
            out[name] = value
    return out


def _coerce(annotation: Any, value: Any) -> Any:
    text = str(annotation)
    if "bool" in text:
        if isinstance(value, bool):
            return value
        return str(value).strip().lower() in {"1", "true", "yes", "on"}
    if "float" in text:
        return float(value)
    if "int" in text and "tuple" not in text:
        return int(value)
    if "tuple" in text:
        if isinstance(value, str):
            return tuple(part.strip() for part in value.split(",") if part.strip())
        return tuple(value)
    return value


_ACTIVE: Settings | None = None


def get_settings(**overrides: Any) -> Settings:
    """Return the process-wide settings object, creating it on first use."""
    global _ACTIVE
    if _ACTIVE is None or overrides:
        _ACTIVE = Settings.load(**overrides)
    return _ACTIVE


def set_settings(settings: Settings) -> None:
    global _ACTIVE
    _ACTIVE = settings


__all__ = [
    "ALPHAFOLD_API_URL",
    "ALPHAFOLD_FILE_URL",
    "DEFAULT_CONFIG_DIR",
    "ESMATLAS_FOLD_URL",
    "RCSB_DOWNLOAD_URL",
    "RCSB_ENTRY_URL",
    "RCSB_POLYMER_URL",
    "RCSB_SEARCH_URL",
    "UNIPROT_FASTA_URL",
    "UNIPROT_JSON_URL",
    "Settings",
    "get_settings",
    "set_settings",
]
