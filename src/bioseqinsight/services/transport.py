"""Pluggable HTTP transports.

The transport is separated from the retry logic on purpose. It is what makes
the fault-tolerance evaluation in ``benchmarks/scripts/run_fault_injection.py``
possible: :class:`FakeTransport` replays a scripted sequence of status codes
and timeouts, so retry, backoff and fallback behaviour can be measured
deterministically and offline instead of waiting for a public service to fail
by chance.

Three implementations are provided:

``RequestsTransport``
    Used when the optional ``requests`` package is installed.
``UrllibTransport``
    Standard-library fallback, so BioSeqInsight retrieves structures on a
    bare Python installation with no third-party HTTP stack at all.
``FakeTransport``
    Test and benchmark double. Never performs I/O.
"""

from __future__ import annotations

import socket
import time
import urllib.error
import urllib.request
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from typing import Protocol

from .. import __version__
from .errors import TransportError

USER_AGENT = f"BioSeqInsight/{__version__} (+https://github.com/niloy-999/BioSeqInsight)"


@dataclass
class Response:
    """Normalised HTTP response."""

    status_code: int
    text: str
    url: str
    headers: dict[str, str] = field(default_factory=dict)
    elapsed_s: float = 0.0

    @property
    def ok(self) -> bool:
        return 200 <= self.status_code < 300


class Transport(Protocol):
    """Minimal transport interface used by :class:`~.http_client.HttpClient`."""

    name: str

    def request(
        self,
        method: str,
        url: str,
        *,
        data: bytes | None = None,
        headers: dict[str, str] | None = None,
        timeout: float = 30.0,
    ) -> Response:  # pragma: no cover - protocol definition
        ...


class UrllibTransport:
    """Standard-library transport. No third-party dependency required."""

    name = "urllib"

    def request(
        self,
        method: str,
        url: str,
        *,
        data: bytes | None = None,
        headers: dict[str, str] | None = None,
        timeout: float = 30.0,
    ) -> Response:
        if not url.lower().startswith("https://"):
            raise TransportError(
                f"Refusing a non-HTTPS request to {url!r}. All bundled structural "
                "endpoints use HTTPS.",
                url=url,
            )
        merged = {"User-Agent": USER_AGENT}
        merged.update(headers or {})
        request = urllib.request.Request(url, data=data, headers=merged, method=method)
        started = time.perf_counter()
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                body = response.read()
                return Response(
                    status_code=response.status,
                    text=body.decode("utf-8", errors="replace"),
                    url=url,
                    headers={k.lower(): v for k, v in response.headers.items()},
                    elapsed_s=time.perf_counter() - started,
                )
        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace") if exc.fp else ""
            return Response(
                status_code=exc.code,
                text=body,
                url=url,
                headers={k.lower(): v for k, v in (exc.headers or {}).items()},
                elapsed_s=time.perf_counter() - started,
            )
        except TimeoutError as exc:
            raise TransportError(f"Timed out after {timeout:.0f}s", url=url, timeout=True) from exc
        except urllib.error.URLError as exc:
            reason = getattr(exc, "reason", exc)
            is_timeout = isinstance(reason, (socket.timeout, TimeoutError))
            raise TransportError(str(reason), url=url, timeout=is_timeout) from exc
        except OSError as exc:  # pragma: no cover - platform specific
            raise TransportError(str(exc), url=url) from exc


class RequestsTransport:
    """Transport backed by the optional ``requests`` package."""

    name = "requests"

    def __init__(self) -> None:
        try:
            import requests
        except ImportError as exc:  # pragma: no cover - exercised only without requests
            raise TransportError(
                "The 'requests' package is not installed. "
                "Install it, or use UrllibTransport."
            ) from exc
        self._requests = requests
        self._session = requests.Session()
        self._session.headers.update({"User-Agent": USER_AGENT})

    def request(
        self,
        method: str,
        url: str,
        *,
        data: bytes | None = None,
        headers: dict[str, str] | None = None,
        timeout: float = 30.0,
    ) -> Response:
        if not url.lower().startswith("https://"):
            raise TransportError(f"Refusing a non-HTTPS request to {url!r}.", url=url)
        started = time.perf_counter()
        try:
            response = self._session.request(
                method, url, data=data, headers=headers, timeout=timeout
            )
        except self._requests.exceptions.Timeout as exc:
            raise TransportError(f"Timed out after {timeout:.0f}s", url=url, timeout=True) from exc
        except self._requests.exceptions.RequestException as exc:
            raise TransportError(str(exc), url=url) from exc
        return Response(
            status_code=response.status_code,
            text=response.text,
            url=url,
            headers={k.lower(): v for k, v in response.headers.items()},
            elapsed_s=time.perf_counter() - started,
        )

    def close(self) -> None:  # pragma: no cover - trivial
        self._session.close()


@dataclass
class ScriptedResponse:
    """One programmed outcome for :class:`FakeTransport`."""

    status_code: int = 200
    body: str = ""
    raise_timeout: bool = False
    raise_error: str | None = None
    delay_s: float = 0.0
    headers: dict[str, str] = field(default_factory=dict)


class FakeTransport:
    """Deterministic transport double used by tests and fault injection.

    ``script`` maps a URL substring to an iterable of
    :class:`ScriptedResponse` objects consumed in order; once exhausted the
    last entry repeats. ``default`` is used for URLs that match nothing.
    """

    name = "fake"

    def __init__(
        self,
        script: dict[str, Iterable[ScriptedResponse]] | None = None,
        default: ScriptedResponse | None = None,
        sleeper: Callable[[float], None] | None = None,
    ) -> None:
        self.script = {key: list(value) for key, value in (script or {}).items()}
        self.default = default or ScriptedResponse(status_code=404, body="not found")
        self.calls: list[tuple[str, str]] = []
        self._cursor: dict[str, int] = {}
        self._sleeper = sleeper or (lambda _seconds: None)

    def _match(self, url: str) -> str | None:
        matches = [key for key in self.script if key in url]
        if not matches:
            return None
        return max(matches, key=len)

    def next_for(self, url: str) -> ScriptedResponse:
        key = self._match(url)
        if key is None:
            return self.default
        queue = self.script[key]
        index = self._cursor.get(key, 0)
        self._cursor[key] = index + 1
        return queue[min(index, len(queue) - 1)]

    def request(
        self,
        method: str,
        url: str,
        *,
        data: bytes | None = None,
        headers: dict[str, str] | None = None,
        timeout: float = 30.0,
    ) -> Response:
        self.calls.append((method, url))
        programmed = self.next_for(url)
        if programmed.delay_s:
            self._sleeper(programmed.delay_s)
        if programmed.raise_timeout:
            raise TransportError(f"Timed out after {timeout:.0f}s", url=url, timeout=True)
        if programmed.raise_error:
            raise TransportError(programmed.raise_error, url=url)
        return Response(
            status_code=programmed.status_code,
            text=programmed.body,
            url=url,
            headers=dict(programmed.headers),
            elapsed_s=programmed.delay_s,
        )


def default_transport() -> Transport:
    """Prefer ``requests`` when available, otherwise the stdlib transport."""
    try:
        return RequestsTransport()
    except TransportError:
        return UrllibTransport()


__all__ = [
    "USER_AGENT",
    "FakeTransport",
    "RequestsTransport",
    "Response",
    "ScriptedResponse",
    "Transport",
    "UrllibTransport",
    "default_transport",
]
