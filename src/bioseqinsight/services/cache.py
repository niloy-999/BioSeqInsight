"""On-disk cache for external responses.

Caching matters here for four reasons that all appear in the evaluation:

1. **Reproducibility.** A cached benchmark run can be replayed months later
   and produce byte-identical results, even if the upstream service has
   changed or disappeared.
2. **Courtesy.** Repeated requests for the same accession do not hit public
   infrastructure repeatedly.
3. **Speed.** Cache hits are reported separately in the performance
   benchmark so the two regimes are never conflated.
4. **Offline inspection.** Every entry stores provenance metadata (query,
   source, URL, timestamp, software version, SHA-256 of the payload), so a
   cached structure can be audited without a network connection.

Entries live under ``<cache_dir>/<namespace>/<hash>.payload`` with a sibling
``.json`` metadata file. Nothing is ever silently overwritten: writing an
entry replaces both files atomically via a temporary file and ``os.replace``.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
import tempfile
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .. import __version__
from .errors import CacheError

logger = logging.getLogger("bioseqinsight.cache")


@dataclass
class CacheEntry:
    """Payload plus provenance metadata."""

    key: str
    namespace: str
    query: str
    source: str
    url: str | None
    created_at: str
    software_version: str
    sha256: str
    size_bytes: int
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def age_s(self) -> float:
        try:
            created = datetime.fromisoformat(self.created_at)
        except ValueError:  # pragma: no cover - defensive
            return float("inf")
        if created.tzinfo is None:
            created = created.replace(tzinfo=timezone.utc)
        return (datetime.now(timezone.utc) - created).total_seconds()


def make_key(namespace: str, query: str, **params: Any) -> str:
    """Stable cache key from a namespace, a query and any request parameters."""
    canonical = json.dumps(
        {"namespace": namespace, "query": query, "params": params},
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:40]


class ResponseCache:
    """Filesystem cache with TTL and an explicit refresh switch."""

    def __init__(
        self,
        directory: str | os.PathLike[str],
        ttl_s: float | None = 30 * 24 * 3600,
        enabled: bool = True,
    ) -> None:
        self.directory = Path(directory)
        self.ttl_s = ttl_s
        self.enabled = enabled
        self.hits = 0
        self.misses = 0

    # -- paths --------------------------------------------------------------

    def _dir(self, namespace: str) -> Path:
        return self.directory / namespace

    def _payload_path(self, namespace: str, key: str) -> Path:
        return self._dir(namespace) / f"{key}.payload"

    def _meta_path(self, namespace: str, key: str) -> Path:
        return self._dir(namespace) / f"{key}.json"

    # -- API ----------------------------------------------------------------

    def get(self, namespace: str, key: str, *, refresh: bool = False) -> tuple[str, CacheEntry] | None:
        """Return ``(payload, entry)`` for a live cache hit, else ``None``."""
        if not self.enabled or refresh:
            self.misses += 1
            return None
        payload_path = self._payload_path(namespace, key)
        meta_path = self._meta_path(namespace, key)
        if not payload_path.is_file() or not meta_path.is_file():
            self.misses += 1
            return None
        try:
            entry = CacheEntry(**json.loads(meta_path.read_text(encoding="utf-8")))
            payload = payload_path.read_text(encoding="utf-8")
        except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
            logger.warning("Discarding unreadable cache entry %s/%s: %s", namespace, key, exc)
            self.delete(namespace, key)
            self.misses += 1
            return None
        if self.ttl_s is not None and entry.age_s > self.ttl_s:
            logger.info("Cache entry %s/%s expired (%.0fs old)", namespace, key, entry.age_s)
            self.misses += 1
            return None
        if hashlib.sha256(payload.encode("utf-8")).hexdigest() != entry.sha256:
            logger.warning("Cache entry %s/%s failed its checksum; discarding", namespace, key)
            self.delete(namespace, key)
            self.misses += 1
            return None
        self.hits += 1
        return payload, entry

    def put(
        self,
        namespace: str,
        key: str,
        payload: str,
        *,
        query: str,
        source: str,
        url: str | None = None,
        extra: dict[str, Any] | None = None,
    ) -> CacheEntry:
        """Store a payload with its provenance metadata."""
        if not self.enabled:
            raise CacheError("Cache is disabled; refusing to write.")
        directory = self._dir(namespace)
        try:
            directory.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise CacheError(f"Cannot create cache directory {directory}: {exc}") from exc

        entry = CacheEntry(
            key=key,
            namespace=namespace,
            query=query,
            source=source,
            url=url,
            created_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
            software_version=__version__,
            sha256=hashlib.sha256(payload.encode("utf-8")).hexdigest(),
            size_bytes=len(payload.encode("utf-8")),
            extra=extra or {},
        )
        try:
            _atomic_write(self._payload_path(namespace, key), payload)
            _atomic_write(
                self._meta_path(namespace, key),
                json.dumps(asdict(entry), indent=2, sort_keys=True),
            )
        except OSError as exc:
            raise CacheError(f"Cannot write cache entry {namespace}/{key}: {exc}") from exc
        return entry

    def delete(self, namespace: str, key: str) -> bool:
        removed = False
        for path in (self._payload_path(namespace, key), self._meta_path(namespace, key)):
            try:
                path.unlink()
                removed = True
            except FileNotFoundError:
                continue
            except OSError as exc:  # pragma: no cover - permissions
                logger.warning("Cannot remove %s: %s", path, exc)
        return removed

    def entries(self, namespace: str | None = None) -> list[CacheEntry]:
        """List metadata for every cached entry, newest first."""
        roots = [self._dir(namespace)] if namespace else _subdirectories(self.directory)
        found: list[CacheEntry] = []
        for root in roots:
            if not root.is_dir():
                continue
            for meta in sorted(root.glob("*.json")):
                try:
                    found.append(CacheEntry(**json.loads(meta.read_text(encoding="utf-8"))))
                except (OSError, TypeError, ValueError):  # pragma: no cover - defensive
                    continue
        return sorted(found, key=lambda e: e.created_at, reverse=True)

    def clear(self, namespace: str | None = None) -> int:
        """Delete cached entries; returns how many payloads were removed."""
        roots = [self._dir(namespace)] if namespace else _subdirectories(self.directory)
        removed = 0
        for root in roots:
            if not root.is_dir():
                continue
            removed += len(list(root.glob("*.payload")))
            shutil.rmtree(root, ignore_errors=True)
        return removed

    def stats(self) -> dict[str, Any]:
        total = self.hits + self.misses
        return {
            "enabled": self.enabled,
            "directory": str(self.directory),
            "ttl_s": self.ttl_s,
            "hits": self.hits,
            "misses": self.misses,
            "hit_rate_percent": (100.0 * self.hits / total) if total else 0.0,
            "entries": len(self.entries()),
        }


def _subdirectories(root: Path) -> list[Path]:
    if not root.is_dir():
        return []
    return [p for p in sorted(root.iterdir()) if p.is_dir()]


def _atomic_write(path: Path, text: str) -> None:
    handle = tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", delete=False, dir=str(path.parent), suffix=".tmp"
    )
    try:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    finally:
        handle.close()
    os.replace(handle.name, path)


class NullCache(ResponseCache):
    """A cache that never stores anything, for ``--no-cache`` runs."""

    def __init__(self) -> None:
        super().__init__(directory=Path(tempfile.gettempdir()) / "bioseqinsight-null", enabled=False)

    def get(self, namespace: str, key: str, *, refresh: bool = False):
        self.misses += 1
        return None

    def put(self, *args, **kwargs) -> CacheEntry:
        return CacheEntry(
            key=kwargs.get("key", ""),
            namespace=args[0] if args else "",
            query=kwargs.get("query", ""),
            source=kwargs.get("source", ""),
            url=kwargs.get("url"),
            created_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
            software_version=__version__,
            sha256="",
            size_bytes=0,
        )


def timed(callable_, *args, **kwargs) -> tuple[Any, float]:
    """Run ``callable_`` and return ``(result, elapsed_seconds)``."""
    started = time.perf_counter()
    result = callable_(*args, **kwargs)
    return result, time.perf_counter() - started


__all__ = ["CacheEntry", "NullCache", "ResponseCache", "make_key", "timed"]
