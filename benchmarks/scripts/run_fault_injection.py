#!/usr/bin/env python3
"""Fault-tolerance benchmark (offline, deterministic, no network required).

The v1.0 manuscript reported that 9 of 25 ESM Atlas requests returned HTTP
504 and stopped there. Reporting an outage is not an engineering result. This
script measures what BioSeqInsight *does* about it.

A scripted transport replays a defined failure pattern for every request, so
the same numbers come out on every machine and in every CI run, without
waiting for a public service to fail by chance and without adding load to it.
Twelve scenarios are covered: clean success, single and repeated transient
failures, each retryable status code, permanent failures, connection errors,
read timeouts, a rate-limited endpoint honouring ``Retry-After``, and the
case that motivated payload validation - a gateway answering HTTP 200 with an
error page instead of coordinates.

For each scenario the script reports how many requests succeeded, how many
succeeded only after an initial failure (recoveries), how many were
unrecoverable, and the mean number of attempts.

Usage
-----
    python benchmarks/scripts/run_fault_injection.py
    python benchmarks/scripts/run_fault_injection.py --repeats 200 --out results.json
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from bioseqinsight import __version__  # noqa: E402
from bioseqinsight.services.http_client import (  # noqa: E402
    HttpClient,
    RetryPolicy,
    summarise_attempts,
)
from bioseqinsight.services.logging_setup import configure_logging  # noqa: E402
from bioseqinsight.services.transport import FakeTransport, ScriptedResponse  # noqa: E402

URL = "https://api.example.org/fold"
COORDINATES = "ATOM      1  N   MET A   1      10.000  10.000  10.000  1.00 90.00           N\nEND\n"

#: name -> (description, scripted responses, expected outcome)
SCENARIOS: dict[str, tuple[str, list[ScriptedResponse], str]] = {
    "healthy_service": (
        "Service answers correctly on the first attempt",
        [ScriptedResponse(200, COORDINATES)],
        "success_first_try",
    ),
    "single_504": (
        "One gateway timeout, then success - the common ESM Atlas pattern",
        [ScriptedResponse(504, "Gateway Time-out"), ScriptedResponse(200, COORDINATES)],
        "recovered",
    ),
    "double_504": (
        "Two gateway timeouts, then success",
        [
            ScriptedResponse(504, "Gateway Time-out"),
            ScriptedResponse(504, "Gateway Time-out"),
            ScriptedResponse(200, COORDINATES),
        ],
        "recovered",
    ),
    "persistent_504": (
        "Service is down for the whole retry budget",
        [ScriptedResponse(504, "Gateway Time-out")],
        "unrecoverable",
    ),
    "http_500": (
        "Internal server error, then success",
        [ScriptedResponse(500, "Internal Server Error"), ScriptedResponse(200, COORDINATES)],
        "recovered",
    ),
    "http_502": (
        "Bad gateway, then success",
        [ScriptedResponse(502, "Bad Gateway"), ScriptedResponse(200, COORDINATES)],
        "recovered",
    ),
    "http_503": (
        "Service unavailable, then success",
        [ScriptedResponse(503, "Service Unavailable"), ScriptedResponse(200, COORDINATES)],
        "recovered",
    ),
    "rate_limited_429": (
        "Rate limited with a Retry-After header, then success",
        [
            ScriptedResponse(429, "Too Many Requests", headers={"retry-after": "1"}),
            ScriptedResponse(200, COORDINATES),
        ],
        "recovered",
    ),
    "not_found_404": (
        "Resource genuinely absent - must not be retried",
        [ScriptedResponse(404, "Not Found")],
        "permanent_no_retry",
    ),
    "forbidden_403": (
        "Permanent client error - must not be retried",
        [ScriptedResponse(403, "Forbidden")],
        "permanent_no_retry",
    ),
    "read_timeout": (
        "Connection established but no response within the timeout",
        [ScriptedResponse(raise_timeout=True), ScriptedResponse(200, COORDINATES)],
        "recovered",
    ),
    "connection_error": (
        "DNS or TCP failure for the whole retry budget",
        [ScriptedResponse(raise_error="Name or service not known")],
        "unrecoverable",
    ),
    "html_error_page_as_200": (
        "Gateway answers HTTP 200 with an HTML error page instead of coordinates",
        [
            ScriptedResponse(200, "<html><body>503 backend unavailable</body></html>"),
            ScriptedResponse(200, COORDINATES),
        ],
        "recovered",
    ),
}


def run_scenario(responses: list[ScriptedResponse], max_attempts: int, repeats: int) -> dict:
    outcomes = []
    for _ in range(repeats):
        transport = FakeTransport({"example.org": list(responses)})
        client = HttpClient(
            transport=transport,
            policy=RetryPolicy(
                max_attempts=max_attempts,
                backoff_base_s=0.0,
                backoff_max_s=0.0,
                jitter=False,
            ),
            timeout_s=5.0,
            sleeper=lambda _seconds: None,  # backoff is measured, not slept through
        )
        outcomes.append(client.post(URL, data=b"MKTAY", expect_text="ATOM"))
    stats = summarise_attempts(outcomes)
    stats["total_http_calls"] = sum(len(o.attempts) for o in outcomes)
    return stats


def classify(stats: dict) -> str:
    if stats["successful"] and stats["first_attempt_success"] == stats["requests"]:
        return "success_first_try"
    if stats["successful"] == stats["requests"] and stats["recovered"]:
        return "recovered"
    if stats["successful"] == 0 and stats["mean_attempts"] == 1.0:
        return "permanent_no_retry"
    return "unrecoverable"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--repeats", type=int, default=100, help="requests per scenario")
    parser.add_argument("--max-attempts", type=int, default=3, help="retry budget")
    parser.add_argument(
        "--out", default="benchmarks/results/fault_injection.json", help="output JSON path"
    )
    args = parser.parse_args(argv)

    # The point of this run is the aggregate table, not a per-attempt log of
    # deliberately induced failures.
    configure_logging(console=False, level="ERROR")

    report = {
        "software_version": __version__,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "repeats_per_scenario": args.repeats,
        "max_attempts": args.max_attempts,
        "note": (
            "Deterministic scripted-transport experiment. No external service was "
            "contacted. Backoff delays are computed by the retry policy but not "
            "slept through, so wall-clock timings are not reported here."
        ),
        "scenarios": {},
    }

    width = max(len(name) for name in SCENARIOS)
    print(
        f"{'scenario':<{width}}  {'success':>8} {'1st try':>8} {'recovered':>10} "
        f"{'failed':>7} {'attempts':>9}  verdict"
    )
    print("-" * (width + 62))

    for name, (description, responses, expected) in SCENARIOS.items():
        stats = run_scenario(responses, args.max_attempts, args.repeats)
        verdict = classify(stats)
        report["scenarios"][name] = {
            "description": description,
            "expected": expected,
            "observed": verdict,
            "as_expected": verdict == expected,
            **stats,
        }
        flag = "ok" if verdict == expected else f"UNEXPECTED (wanted {expected})"
        print(
            f"{name:<{width}}  {stats['successful']:>8} {stats['first_attempt_success']:>8} "
            f"{stats['recovered']:>10} {stats['unrecoverable']:>7} "
            f"{stats['mean_attempts']:>9.2f}  {flag}"
        )

    transient = [
        name
        for name, (_, _, expected) in SCENARIOS.items()
        if expected in ("recovered", "unrecoverable")
    ]
    recovered = sum(report["scenarios"][n]["recovered"] for n in transient)
    initial_failures = sum(report["scenarios"][n]["initial_failures"] for n in transient)
    report["aggregate"] = {
        "transient_scenarios": len(transient),
        "requests_with_an_initial_failure": initial_failures,
        "recovered": recovered,
        "recovery_rate_percent": (
            round(100.0 * recovered / initial_failures, 1) if initial_failures else 0.0
        ),
        "all_scenarios_behaved_as_specified": all(
            entry["as_expected"] for entry in report["scenarios"].values()
        ),
    }

    print("\nAggregate over scenarios with transient or permanent faults:")
    print(json.dumps(report["aggregate"], indent=2))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nWritten to {out}", file=sys.stderr)
    return 0 if report["aggregate"]["all_scenarios_behaved_as_specified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
