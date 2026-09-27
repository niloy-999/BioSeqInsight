"""Exception hierarchy for the service layer.

Every failure mode that BioSeqInsight can encounter when talking to an
external resource maps onto exactly one of these classes. The benchmark
harness classifies outcomes by exception type, so this hierarchy is part of
the measurement instrument, not just error plumbing.
"""

from __future__ import annotations


class BioSeqInsightError(Exception):
    """Base class for all BioSeqInsight errors."""


class ConfigurationError(BioSeqInsightError):
    """Invalid or contradictory configuration."""


class TransportError(BioSeqInsightError):
    """Network-level failure: DNS, TCP, TLS or a read timeout."""

    def __init__(self, message: str, *, url: str | None = None, timeout: bool = False):
        super().__init__(message)
        self.url = url
        self.timeout = timeout


class HTTPStatusError(BioSeqInsightError):
    """The server answered, but with a status code we cannot use."""

    def __init__(self, status_code: int, message: str = "", *, url: str | None = None, body: str = ""):
        super().__init__(message or f"HTTP {status_code}")
        self.status_code = status_code
        self.url = url
        self.body = body

    @property
    def retryable(self) -> bool:
        return self.status_code in RETRYABLE_STATUS


class ServiceUnavailable(BioSeqInsightError):
    """A provider exhausted its retry budget without succeeding."""

    def __init__(self, service: str, attempts: int, last_error: str = ""):
        super().__init__(
            f"{service} did not return a usable response after {attempts} attempt(s). "
            f"Last error: {last_error or 'unknown'}"
        )
        self.service = service
        self.attempts = attempts
        self.last_error = last_error


class NotFoundError(BioSeqInsightError):
    """The resource genuinely does not hold a record for this query.

    Distinct from :class:`ServiceUnavailable`: a 404 from AlphaFold DB means
    "no model for this accession", which is a *result*, whereas a 504 means
    "ask again later", which is not.
    """


class ValidationFailure(BioSeqInsightError):
    """A response was received but failed structural or content validation."""


class CacheError(BioSeqInsightError):
    """The on-disk cache could not be read or written."""


RETRYABLE_STATUS = frozenset({408, 409, 425, 429, 500, 502, 503, 504})
PERMANENT_STATUS = frozenset({400, 401, 403, 404, 405, 410, 422})


__all__ = [
    "PERMANENT_STATUS",
    "RETRYABLE_STATUS",
    "BioSeqInsightError",
    "CacheError",
    "ConfigurationError",
    "HTTPStatusError",
    "NotFoundError",
    "ServiceUnavailable",
    "TransportError",
    "ValidationFailure",
]
