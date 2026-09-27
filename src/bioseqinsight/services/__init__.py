"""Cross-cutting infrastructure: transport, retry, cache and logging."""

from .cache import NullCache, ResponseCache
from .errors import (
    BioSeqInsightError,
    HTTPStatusError,
    NotFoundError,
    ServiceUnavailable,
    TransportError,
    ValidationFailure,
)
from .http_client import HttpClient, RequestOutcome, RetryPolicy, summarise_attempts
from .logging_setup import configure_logging, get_logger
from .transport import FakeTransport, Response, ScriptedResponse, UrllibTransport

__all__ = [
    "BioSeqInsightError",
    "FakeTransport",
    "HTTPStatusError",
    "HttpClient",
    "NotFoundError",
    "NullCache",
    "RequestOutcome",
    "Response",
    "ResponseCache",
    "RetryPolicy",
    "ScriptedResponse",
    "ServiceUnavailable",
    "TransportError",
    "UrllibTransport",
    "ValidationFailure",
    "configure_logging",
    "get_logger",
    "summarise_attempts",
]
