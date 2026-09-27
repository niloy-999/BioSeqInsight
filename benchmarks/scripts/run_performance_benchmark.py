#!/usr/bin/env python3
"""Performance and scalability benchmark (offline, no network required).

Measures the parts of BioSeqInsight whose speed is under its own control:

* local nucleotide and protein analysis, over batches of increasing size;
* the pairwise alignment that underpins structure validation, as a function
  of sequence length;
* cache write and cache read latency, and the speed-up a cache hit provides
  over a (simulated) 2-second remote call;
* export throughput for CSV, JSON and HTML;
* application start-up cost, measured as the time to import the package and
  construct the analysis stack.

Remote response times are deliberately **not** measured here. They are a
property of RCSB, EBI and ESM Atlas rather than of this software, they vary
with time of day and location, and mixing them into a scalability curve would
make the curve meaningless. Live service timings are reported separately by
``run_structure_benchmark.py``.

Every input is generated from a fixed random seed, so the workload is
identical on every machine and the numbers can be compared across platforms.

Usage
-----
    python benchmarks/scripts/run_performance_benchmark.py
    python benchmarks/scripts/run_performance_benchmark.py --sizes 10,100,1000 --repeats 5
"""

from __future__ import annotations

import argparse
import json
import platform
import random
import statistics
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from bioseqinsight import __version__  # noqa: E402
from bioseqinsight.config.settings import Settings  # noqa: E402
from bioseqinsight.core.alignment import semi_global  # noqa: E402
from bioseqinsight.core.alphabet import FastaRecord  # noqa: E402
from bioseqinsight.core.protein import ProteinAnalyzer  # noqa: E402
from bioseqinsight.core.sequence import SequenceAnalyzer  # noqa: E402
from bioseqinsight.io.export import build_html_report, rows_to_delimited, to_json  # noqa: E402
from bioseqinsight.services.cache import ResponseCache, make_key  # noqa: E402
from bioseqinsight.workflows.batch import BatchOptions, BatchRunner  # noqa: E402

SEED = 20260921
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
SIMULATED_REMOTE_CALL_S = 2.0  # median ESM Atlas response time reported in v1.0


def make_dna(rng: random.Random, length: int) -> str:
    return "ATG" + "".join(rng.choice("ACGT") for _ in range(length - 6)) + "TAA"


def make_protein(rng: random.Random, length: int) -> str:
    return "M" + "".join(rng.choice(AMINO_ACIDS) for _ in range(length - 1))


def timed(function, *args, repeats: int = 1, **kwargs) -> tuple[float, float]:
    """Run ``function`` ``repeats`` times; return (median seconds, stdev)."""
    samples = []
    for _ in range(repeats):
        started = time.perf_counter()
        function(*args, **kwargs)
        samples.append(time.perf_counter() - started)
    return statistics.median(samples), (statistics.stdev(samples) if len(samples) > 1 else 0.0)


def benchmark_batch(sizes: list[int], repeats: int) -> list[dict]:
    rng = random.Random(SEED)
    rows = []
    with TemporaryDirectory() as tmp:
        settings = Settings.from_dict(
            {
                "download_dir": f"{tmp}/s",
                "cache_dir": f"{tmp}/c",
                "projects_dir": f"{tmp}/p",
                "log_file": f"{tmp}/l/x.log",
                "offline": True,
            }
        )
        runner = BatchRunner(settings)
        for size in sizes:
            dna = [
                FastaRecord(f"dna{i}", "", make_dna(rng, 900)) for i in range(size // 2)
            ]
            protein = [
                FastaRecord(f"prot{i}", "", make_protein(rng, 300))
                for i in range(size - size // 2)
            ]
            records = dna + protein
            median, stdev = timed(
                runner.run, records, BatchOptions(include_structures=False), repeats=repeats
            )
            rows.append(
                {
                    "sequences": size,
                    "median_s": round(median, 4),
                    "stdev_s": round(stdev, 4),
                    "per_sequence_ms": round(1000.0 * median / size, 3),
                    "sequences_per_second": round(size / median, 1),
                }
            )
            print(
                f"  batch {size:>5} sequences: {median:>7.3f} s "
                f"({1000.0 * median / size:.2f} ms each, {size / median:.0f}/s)",
                file=sys.stderr,
            )
    return rows


def benchmark_single(repeats: int) -> dict:
    rng = random.Random(SEED)
    results = {}
    for length in (300, 3000, 30000):
        sequence = make_dna(rng, length)
        median, _ = timed(
            lambda s=sequence: SequenceAnalyzer(s).analyze(min_orf_aa=30), repeats=repeats
        )
        results[f"dna_{length}nt_ms"] = round(1000 * median, 3)
    for length in (100, 300, 1000):
        sequence = make_protein(rng, length)
        median, _ = timed(lambda s=sequence: ProteinAnalyzer(s).analyze(), repeats=repeats)
        results[f"protein_{length}aa_ms"] = round(1000 * median, 3)
    return results


def benchmark_alignment(repeats: int) -> list[dict]:
    rng = random.Random(SEED)
    rows = []
    for length in (76, 150, 300, 600):
        query = make_protein(rng, length)
        subject = query[: int(length * 0.9)]
        median, _ = timed(lambda q=query, s=subject: semi_global(q, s), repeats=repeats)
        rows.append(
            {
                "query_length": length,
                "subject_length": len(subject),
                "median_ms": round(1000 * median, 2),
                "cells": (length + 1) * (len(subject) + 1),
            }
        )
        print(f"  alignment {length:>4} x {len(subject):>4}: {1000 * median:>7.1f} ms", file=sys.stderr)
    return rows


def benchmark_cache(repeats: int) -> dict:
    payload = "ATOM  " + "x" * 200_000  # roughly the size of a real coordinate file
    with TemporaryDirectory() as tmp:
        cache = ResponseCache(Path(tmp) / "cache", ttl_s=None)
        key = make_key("bench", "P24941")
        write_median, _ = timed(
            lambda: cache.put("bench", key, payload, query="P24941", source="bench"),
            repeats=repeats,
        )
        read_median, _ = timed(lambda: cache.get("bench", key), repeats=repeats)
    return {
        "payload_bytes": len(payload),
        "write_ms": round(1000 * write_median, 3),
        "read_ms": round(1000 * read_median, 3),
        "simulated_remote_call_s": SIMULATED_REMOTE_CALL_S,
        "speedup_versus_remote_call": round(SIMULATED_REMOTE_CALL_S / read_median, 0),
        "note": (
            "The remote call is a constant taken from the v1.0 benchmark's median "
            "ESM Atlas response time, not a measurement made here. The point of the "
            "ratio is the order of magnitude: a cache hit is local disk I/O."
        ),
    }


def benchmark_export(repeats: int) -> dict:
    rng = random.Random(SEED)
    with TemporaryDirectory() as tmp:
        settings = Settings.from_dict(
            {"download_dir": f"{tmp}/s", "cache_dir": f"{tmp}/c", "log_file": f"{tmp}/l/x.log",
             "offline": True}
        )
        records = [FastaRecord(f"p{i}", "", make_protein(rng, 300)) for i in range(200)]
        rows, _ = BatchRunner(settings).run(records)
        csv_median, _ = timed(lambda: rows_to_delimited(rows), repeats=repeats)
        json_median, _ = timed(lambda: to_json([r.to_dict() for r in rows]), repeats=repeats)
        html_median, _ = timed(lambda: build_html_report(rows), repeats=repeats)
    return {
        "rows": len(rows),
        "csv_ms": round(1000 * csv_median, 2),
        "json_ms": round(1000 * json_median, 2),
        "html_ms": round(1000 * html_median, 2),
    }


def benchmark_startup() -> dict:
    """Cost of standing up the analysis stack in a fresh interpreter."""
    import subprocess

    code = (
        "import sys, time; sys.path.insert(0, r'%s');"
        "t=time.perf_counter();"
        "import bioseqinsight;"
        "from bioseqinsight.core.sequence import SequenceAnalyzer;"
        "from bioseqinsight.core.protein import ProteinAnalyzer;"
        "SequenceAnalyzer('ATGGCTTAA').analyze();"
        "print(round((time.perf_counter()-t)*1000, 1))" % str(ROOT / "src")
    )
    samples = []
    for _ in range(5):
        output = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, check=True
        )
        samples.append(float(output.stdout.strip()))
    return {
        "import_and_first_analysis_ms": round(statistics.median(samples), 1),
        "samples_ms": samples,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sizes", default="10,100,500,1000")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--out", default="benchmarks/results/performance.json")
    args = parser.parse_args(argv)
    sizes = [int(value) for value in args.sizes.split(",") if value.strip()]

    print("Single-sequence analysis...", file=sys.stderr)
    single = benchmark_single(args.repeats)
    print("Batch scalability...", file=sys.stderr)
    batch = benchmark_batch(sizes, args.repeats)
    print("Alignment (structure validation)...", file=sys.stderr)
    alignment = benchmark_alignment(args.repeats)
    print("Cache...", file=sys.stderr)
    cache = benchmark_cache(max(args.repeats, 5))
    print("Export...", file=sys.stderr)
    export = benchmark_export(args.repeats)
    print("Start-up...", file=sys.stderr)
    startup = benchmark_startup()

    # Linearity check: if the per-sequence cost is roughly flat across the
    # range, throughput scales linearly with batch size.
    per_sequence = [row["per_sequence_ms"] for row in batch]
    linearity = (max(per_sequence) / min(per_sequence)) if min(per_sequence) > 0 else None

    report = {
        "software_version": __version__,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "environment": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "processor": platform.processor() or "unknown",
        },
        "seed": SEED,
        "repeats": args.repeats,
        "single_sequence": single,
        "batch_scalability": batch,
        "batch_per_sequence_cost_ratio": round(linearity, 2) if linearity else None,
        "alignment": alignment,
        "cache": cache,
        "export": export,
        "startup": startup,
        "note": (
            "Local computation only. External service latency is a property of the "
            "remote resources and is reported by run_structure_benchmark.py."
        ),
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print("\n== Summary ==")
    print(json.dumps(
        {
            "single_sequence": single,
            "batch_scalability": batch,
            "batch_per_sequence_cost_ratio": report["batch_per_sequence_cost_ratio"],
            "cache_read_ms": cache["read_ms"],
            "startup_ms": startup["import_and_first_analysis_ms"],
        },
        indent=2,
    ))
    print(f"\nWritten to {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
