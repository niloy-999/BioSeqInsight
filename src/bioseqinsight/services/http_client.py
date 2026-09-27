"""Fault-tolerant HTTP client.

The v1.0 manuscript reported that 9 of 25 ESM Atlas requests returned HTTP
504. v2.0 treats that as an engineering requirement rather than an
observation: every external call goes through this client, which applies a
bounded retry policy with exponential backoff and full jitter, honours
``Retry-After``, distinguishes permanent from transient failures, and records
each individual attempt so that recovery can be measured.

Every attempt is returned to the caller as an
:class:`~bioseqinsight.models.results.AttemptRecord`, which is what lets the
benchmark report *initial failures*, *recovered failures*, *unrecoverable
failures* and *mean recovery time* rather than a single success count.
"""

from __future__ import annotations

import logging
import random
import time
from dataclasses import dataclass, field

from ..models.results import AttemptRecord
from .errors import (
    PERMANENT_STATUS,
    RETRYABLE_STATUS,
    HTTPStatusError,
    NotFoundError,
    ServiceUnavailable,
    TransportError,
)
from .transport import Response, Transport, default_transport

logger = logging.getLogger("bioseqinsight.http")


@dataclass(frozen=True)
class RetryPolicy:
    """Bounded exponential backoff with full jitter.

    ``max_attempts`` counts the first try, so ``max_attempts=3`` means one
    request plus at most two retries. ``jitter`` randomises the delay across
    ``[0, computed_delay]``; with ``jitter=False`` the delay is deterministic,
    which the tests rely on.
    """

    max_attempts: int = 3
    backoff_base_s: float = 1.0
    backoff_factor: float = 2.0
    backoff_max_s: float = 20.0
    jitter: bool = True
    retry_statuses: frozenset[int] = RETRYABLE_STATUS
    respect_retry_after: bool = True
    max_retry_after_s: float = 30.0

    def delay_for(self, attempt: int, retry_after: float | None = None) -> float:
        """Seconds to wait before attempt number ``attempt + 1``."""
        if self.respect_retry_after and retry_after is not None:
            return min(max(retry_after, 0.0), self.max_retry_after_s)
        raw = min(self.backoff_base_s * (self.backoff_factor ** (attempt - 1)), self.backoff_max_s)
        if not self.jitter:
            return raw
        return random.uniform(0.0, raw)

    def should_retry(self, status_code: int) -> bool:
        return status_code in self.retry_statuses


@dataclass
class RequestOutcome:
    """A completed request, successful or not, with its full attempt history."""

    response: Response | None
    attempts: list[AttemptRecord] = field(default_factory=list)
    error: Exception | None = None
    total_elapsed_s: float = 0.0

    @property
    def ok(self) -> bool:
        return self.response is not None and self.response.ok

    @property
    def recovered(self) -> bool:
        """Succeeded, but only after at least one failed attempt."""
        return self.ok and len(self.attempts) > 1

    @property
    def status_code(self) -> int | None:
        return self.response.status_code if self.response else None


def _parse_retry_after(response: Response) -> float | None:
    value = response.headers.get("retry-after")
    if not value:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None  # HTTP-date form; fall back to exponential backoff


class HttpClient:
    """Retrying HTTP client shared by every structural provider."""

    def __init__(
        self,
        transport: Transport | None = None,
        policy: RetryPolicy | None = None,
        timeout_s: float = 30.0,
        sleeper=None,
    ) -> None:
        self.transport = transport or default_transport()
        self.policy = policy or RetryPolicy()
        self.timeout_s = timeout_s
        self._sleep = sleeper or time.sleep

    def request(
        self,
        method: str,
        url: str,
        *,
        data: bytes | None = None,
        headers: dict[str, str] | None = None,
        timeout_s: float | None = None,
        expect_text: str | None = None,
    ) -> RequestOutcome:
        """Perform a request with retries and return its full history.

        ``expect_text`` validates the payload: if the substring is absent the
        response is treated as a failed attempt. Structural endpoints
        occasionally answer 200 with an HTML error page, and accepting that
        as a structure is exactly the kind of silent corruption v2.0 exists
        to prevent.
        """
        outcome = RequestOutcome(response=None)
        started = time.perf_counter()
        effective_timeout = timeout_s if timeout_s is not None else self.timeout_s

        for attempt in range(1, self.policy.max_attempts + 1):
            attempt_started = time.perf_counter()
            try:
                response = self.transport.request(
                    method, url, data=data, headers=headers, timeout=effective_timeout
                )
            except TransportError as exc:
                elapsed = time.perf_counter() - attempt_started
                retryable = attempt < self.policy.max_attempts
                outcome.attempts.append(
                    AttemptRecord(
                        attempt=attempt,
                        url=url,
                        status_code=None,
                        error=("timeout: " if exc.timeout else "") + str(exc),
                        elapsed_s=elapsed,
                        retried=retryable,
                    )
                )
                outcome.error = exc
                logger.warning("%s %s attempt %d failed: %s", method, url, attempt, exc)
                if not retryable:
                    break
                self._sleep(self.policy.delay_for(attempt))
                continue

            elapsed = time.perf_counter() - attempt_started
            payload_bad = (
                response.ok and expect_text is not None and expect_text not in response.text
            )
            success = response.ok and not payload_bad

            error_text: str | None = None
            if payload_bad:
                error_text = (
                    f"HTTP 200 but the payload did not contain {expect_text!r}; "
                    "treating as an invalid response."
                )
            elif not response.ok:
                snippet = response.text[:160].replace("\n", " ").strip()
                error_text = f"HTTP {response.status_code}" + (f" - {snippet}" if snippet else "")

            retryable = (
                not success
                and attempt < self.policy.max_attempts
                and (payload_bad or self.policy.should_retry(response.status_code))
            )
            outcome.attempts.append(
                AttemptRecord(
                    attempt=attempt,
                    url=url,
                    status_code=response.status_code,
                    error=error_text,
                    elapsed_s=elapsed,
                    retried=retryable,
                )
            )

            if success:
                outcome.response = response
                outcome.error = None
                outcome.total_elapsed_s = time.perf_counter() - started
                if attempt > 1:
                    logger.info("%s %s recovered on attempt %d", method, url, attempt)
                return outcome

            outcome.response = response
            outcome.error = HTTPStatusError(
                response.status_code, error_text or "", url=url, body=response.text[:500]
            )
            if response.status_code in PERMANENT_STATUS:
                logger.info("%s %s permanent failure %s", method, url, response.status_code)
                break
            if not retryable:
                break
            self._sleep(self.policy.delay_for(attempt, _parse_retry_after(response)))

        outcome.total_elapsed_s = time.perf_counter() - started
        return outcome

    def get(self, url: str, **kwargs) -> RequestOutcome:
        return self.request("GET", url, **kwargs)

    def post(self, url: str, data: bytes, **kwargs) -> RequestOutcome:
        return self.request("POST", url, data=data, **kwargs)

    def get_text(self, url: str, service: str = "service", **kwargs) -> str:
        """Convenience wrapper that raises instead of returning an outcome."""
        outcome = self.get(url, **kwargs)
        if outcome.ok and outcome.response is not None:
            return outcome.response.text
        if outcome.status_code == 404:
            raise NotFoundError(f"{service} has no record for {url}")
        raise ServiceUnavailable(service, len(outcome.attempts), str(outcome.error or ""))


def summarise_attempts(outcomes: list[RequestOutcome]) -> dict[str, float]:
    """Aggregate statistics used by the reliability benchmark."""
    total = len(outcomes)
    if not total:
        return {
            "requests": 0,
            "successful": 0,
            "first_attempt_success": 0,
            "recovered": 0,
            "unrecoverable": 0,
            "success_rate": 0.0,
            "recovery_rate": 0.0,
            "mean_attempts": 0.0,
            "mean_recovery_time_s": 0.0,
        }
    successful = [o for o in outcomes if o.ok]
    recovered = [o for o in successful if o.recovered]
    first_try = [o for o in successful if len(o.attempts) == 1]
    initial_failures = [o for o in outcomes if o.attempts and o.attempts[0].error]
    recovery_times = [o.total_elapsed_s for o in recovered]
    return {
        "requests": total,
        "successful": len(successful),
        "first_attempt_success": len(first_try),
        "initial_failures": len(initial_failures),
        "recovered": len(recovered),
        "unrecoverable": total - len(successful),
        "success_rate": 100.0 * len(successful) / total,
        "recovery_rate": (
            100.0 * len(recovered) / len(initial_failures) if initial_failures else 0.0
        ),
        "mean_attempts": sum(len(o.attempts) for o in outcomes) / total,
        "mean_recovery_time_s": (
            sum(recovery_times) / len(recovery_times) if recovery_times else 0.0
        ),
    }


__all__ = [
    "HttpClient",
    "RequestOutcome",
    "RetryPolicy",
    "summarise_attempts",
]
