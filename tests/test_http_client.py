"""Retry policy, backoff and attempt recording.

These tests are the offline half of the fault-tolerance evaluation: they
prove the recovery machinery works without waiting for a public service to
fail by chance.
"""

from __future__ import annotations

import pytest

from bioseqinsight.services.errors import TransportError
from bioseqinsight.services.http_client import (
    RetryPolicy,
    summarise_attempts,
)
from bioseqinsight.services.transport import FakeTransport, ScriptedResponse
from conftest import make_client

URL = "https://example.org/resource"


def transport(*responses: ScriptedResponse) -> FakeTransport:
    return FakeTransport({"example.org": list(responses)})


class TestRetryPolicy:
    def test_backoff_grows_exponentially(self):
        policy = RetryPolicy(backoff_base_s=1.0, backoff_factor=2.0, jitter=False)
        assert [policy.delay_for(n) for n in (1, 2, 3)] == [1.0, 2.0, 4.0]

    def test_backoff_is_capped(self):
        policy = RetryPolicy(backoff_base_s=1.0, backoff_max_s=3.0, jitter=False)
        assert policy.delay_for(10) == 3.0

    def test_jitter_stays_within_the_bound(self):
        policy = RetryPolicy(backoff_base_s=4.0, jitter=True)
        assert all(0.0 <= policy.delay_for(1) <= 4.0 for _ in range(50))

    def test_retry_after_header_wins(self):
        policy = RetryPolicy(backoff_base_s=10.0, jitter=False)
        assert policy.delay_for(1, retry_after=2.0) == 2.0

    def test_retry_after_is_capped(self):
        policy = RetryPolicy(max_retry_after_s=5.0, jitter=False)
        assert policy.delay_for(1, retry_after=9999.0) == 5.0

    def test_transient_statuses_are_retryable(self):
        policy = RetryPolicy()
        assert all(policy.should_retry(code) for code in (429, 500, 502, 503, 504))

    def test_permanent_statuses_are_not_retryable(self):
        policy = RetryPolicy()
        assert not any(policy.should_retry(code) for code in (400, 404, 403))


class TestRequestBehaviour:
    def test_success_on_the_first_attempt(self):
        client = make_client(transport(ScriptedResponse(200, "ATOM  ok")))
        outcome = client.get(URL)
        assert outcome.ok
        assert len(outcome.attempts) == 1
        assert not outcome.recovered

    def test_recovery_after_transient_failures(self):
        client = make_client(
            transport(
                ScriptedResponse(504, "gateway timeout"),
                ScriptedResponse(503, "unavailable"),
                ScriptedResponse(200, "ATOM  ok"),
            )
        )
        outcome = client.get(URL)
        assert outcome.ok
        assert outcome.recovered
        assert len(outcome.attempts) == 3
        assert outcome.attempts[0].status_code == 504
        assert outcome.attempts[0].retried is True

    def test_retry_budget_is_respected(self):
        client = make_client(transport(ScriptedResponse(504, "gateway timeout")), attempts=3)
        outcome = client.get(URL)
        assert not outcome.ok
        assert len(outcome.attempts) == 3

    def test_permanent_failure_is_not_retried(self):
        client = make_client(transport(ScriptedResponse(404, "missing")))
        outcome = client.get(URL)
        assert not outcome.ok
        assert len(outcome.attempts) == 1
        assert outcome.status_code == 404

    def test_network_timeout_is_recorded_as_such(self):
        client = make_client(transport(ScriptedResponse(raise_timeout=True)), attempts=2)
        outcome = client.get(URL)
        assert not outcome.ok
        assert all(a.error and a.error.startswith("timeout") for a in outcome.attempts)

    def test_timeout_then_success(self):
        client = make_client(
            transport(ScriptedResponse(raise_timeout=True), ScriptedResponse(200, "ATOM"))
        )
        outcome = client.get(URL)
        assert outcome.ok and outcome.recovered

    def test_connection_error_is_recorded(self):
        client = make_client(transport(ScriptedResponse(raise_error="DNS failure")), attempts=1)
        outcome = client.get(URL)
        assert not outcome.ok
        assert "DNS failure" in outcome.attempts[0].error

    def test_payload_validation_rejects_a_200_error_page(self):
        # A gateway that answers 200 with an HTML error page must not be
        # accepted as a structure.
        client = make_client(
            transport(
                ScriptedResponse(200, "<html>503 Service Unavailable</html>"),
                ScriptedResponse(200, "ATOM      1  N   MET A   1"),
            )
        )
        outcome = client.get(URL, expect_text="ATOM")
        assert outcome.ok
        assert outcome.recovered
        assert "did not contain" in outcome.attempts[0].error

    def test_post_sends_a_body(self):
        fake = transport(ScriptedResponse(200, "ATOM"))
        make_client(fake).post(URL, data=b"MKTAY")
        assert fake.calls == [("POST", URL)]

    def test_non_https_is_refused_by_the_stdlib_transport(self):
        from bioseqinsight.services.transport import UrllibTransport

        with pytest.raises(TransportError, match="non-HTTPS"):
            UrllibTransport().request("GET", "http://example.org")


class TestSummaries:
    def test_empty_input(self):
        assert summarise_attempts([])["requests"] == 0

    def test_mixed_outcomes_are_counted(self):
        first = make_client(transport(ScriptedResponse(200, "ATOM"))).get(URL)
        recovered = make_client(
            transport(ScriptedResponse(504, "x"), ScriptedResponse(200, "ATOM"))
        ).get(URL)
        failed = make_client(transport(ScriptedResponse(504, "x")), attempts=2).get(URL)

        stats = summarise_attempts([first, recovered, failed])
        assert stats["requests"] == 3
        assert stats["successful"] == 2
        assert stats["first_attempt_success"] == 1
        assert stats["recovered"] == 1
        assert stats["unrecoverable"] == 1
        assert stats["initial_failures"] == 2
        assert stats["recovery_rate"] == pytest.approx(50.0)
        assert stats["success_rate"] == pytest.approx(200 / 3)


class TestFakeTransport:
    def test_last_scripted_response_repeats_once_exhausted(self):
        fake = transport(ScriptedResponse(500, "a"), ScriptedResponse(404, "b"))
        client = make_client(fake, attempts=1)
        assert client.get(URL).status_code == 500
        assert client.get(URL).status_code == 404
        assert client.get(URL).status_code == 404

    def test_unmatched_urls_get_the_default(self):
        fake = FakeTransport({}, default=ScriptedResponse(418, "teapot"))
        assert make_client(fake, attempts=1).get(URL).status_code == 418

    def test_longest_matching_key_wins(self):
        fake = FakeTransport(
            {
                "example.org": [ScriptedResponse(200, "generic")],
                "example.org/specific": [ScriptedResponse(201, "specific")],
            }
        )
        client = make_client(fake, attempts=1)
        assert client.get("https://example.org/specific").status_code == 201
