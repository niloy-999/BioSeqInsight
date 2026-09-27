"""Logging configuration.

Research software needs an audit trail: what was requested, from which
resource, how long it took and what came back. The log is also the raw
material for the reliability benchmark, so its format is stable and
machine-parseable.

Privacy note: full sequences are never written to the log. Sequence inputs
are recorded as a length plus a short SHA-256 prefix, which is enough to
correlate runs without persisting potentially unpublished sequence data to
disk.
"""

from __future__ import annotations

import hashlib
import json
import logging
import logging.handlers
from pathlib import Path

LOGGER_NAME = "bioseqinsight"
_CONFIGURED = False

TEXT_FORMAT = "%(asctime)s %(levelname)-7s %(name)s: %(message)s"


class JsonFormatter(logging.Formatter):
    """One JSON object per line, for machine-readable audit logs."""

    def format(self, record: logging.LogRecord) -> str:
        payload = {
            "time": self.formatTime(record, "%Y-%m-%dT%H:%M:%S"),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
        }
        for key in ("operation", "resource", "query_ref", "status", "elapsed_s", "version"):
            value = getattr(record, key, None)
            if value is not None:
                payload[key] = value
        if record.exc_info:
            payload["exception"] = self.formatException(record.exc_info)
        return json.dumps(payload, ensure_ascii=False)


def sequence_reference(sequence: str) -> str:
    """Non-reversible reference for a sequence: ``len=298:sha256=1f4a9c22``."""
    digest = hashlib.sha256((sequence or "").encode("utf-8")).hexdigest()[:8]
    return f"len={len(sequence or '')}:sha256={digest}"


def configure_logging(
    log_file: str | Path | None = None,
    level: str = "INFO",
    json_format: bool = False,
    console: bool = True,
    max_bytes: int = 2_000_000,
    backup_count: int = 3,
) -> logging.Logger:
    """Configure the ``bioseqinsight`` logger. Safe to call more than once."""
    global _CONFIGURED
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(getattr(logging, str(level).upper(), logging.INFO))
    logger.propagate = False

    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()

    formatter: logging.Formatter = JsonFormatter() if json_format else logging.Formatter(TEXT_FORMAT)

    if console:
        stream = logging.StreamHandler()
        stream.setFormatter(formatter)
        logger.addHandler(stream)

    if log_file:
        path = Path(log_file)
        path.parent.mkdir(parents=True, exist_ok=True)
        rotating = logging.handlers.RotatingFileHandler(
            path, maxBytes=max_bytes, backupCount=backup_count, encoding="utf-8"
        )
        rotating.setFormatter(JsonFormatter() if json_format else logging.Formatter(TEXT_FORMAT))
        logger.addHandler(rotating)

    if not logger.handlers:
        logger.addHandler(logging.NullHandler())

    _CONFIGURED = True
    return logger


def get_logger(name: str = "") -> logging.Logger:
    """Return a child of the package logger, configuring a default if needed."""
    if not _CONFIGURED:
        configure_logging(console=False)
    return logging.getLogger(f"{LOGGER_NAME}.{name}" if name else LOGGER_NAME)


def log_operation(
    logger: logging.Logger,
    operation: str,
    *,
    resource: str = "local",
    query_ref: str = "",
    status: str = "ok",
    elapsed_s: float | None = None,
    level: int = logging.INFO,
) -> None:
    """Emit a structured operation record."""
    logger.log(
        level,
        "%s on %s: %s",
        operation,
        resource,
        status,
        extra={
            "operation": operation,
            "resource": resource,
            "query_ref": query_ref,
            "status": status,
            "elapsed_s": round(elapsed_s, 4) if elapsed_s is not None else None,
        },
    )


__all__ = [
    "LOGGER_NAME",
    "JsonFormatter",
    "configure_logging",
    "get_logger",
    "log_operation",
    "sequence_reference",
]
