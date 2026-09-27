#!/usr/bin/env python3
"""Structure-resource benchmark (requires network access).

What this script does
---------------------
1. Reads the curated candidate list in ``benchmarks/data/benchmark_candidates.csv``.
2. Resolves every accession against UniProt, recording the canonical sequence
   and its true length. Accessions that do not resolve are excluded and
   listed, so the resolved set is always reproducible and auditable.
3. Assigns each resolved protein to a length stratum (short < 150 aa,
   medium 150-300 aa, long > 300 aa) **from the fetched length**, never from a
   figure typed into the dataset by hand.
4. Queries every configured structural resource for every protein, validates
   each returned structure against the canonical sequence, and records the
   mapping level, identity, coverage, attempt count and response time.
5. Writes a per-protein table, a machine-readable JSON record and a summary.

Why the strata are computed rather than stored
----------------------------------------------
A benchmark table that states a sequence length is only as trustworthy as
whoever typed it. Fetching the length at run time means the published
stratification cannot silently disagree with the sequences actually analysed,
and any drift in a UniProt entry between runs shows up as a changed resolved
file rather than as an invisible error.

Usage
-----
    python benchmarks/scripts/run_structure_benchmark.py \
        --out benchmarks/results/structure_benchmark_2026-09-21

    # a quick smoke run over the first five proteins
    python benchmarks/scripts/run_structure_benchmark.py --limit 5
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from bioseqinsight import __version__  # noqa: E402
from bioseqinsight.config.settings import Settings  # noqa: E402
from bioseqinsight.models.results import MappingLevel  # noqa: E402
from bioseqinsight.services.errors import BioSeqInsightError  # noqa: E402
from bioseqinsight.services.logging_setup import configure_logging  # noqa: E402
from bioseqinsight.structures.manager import StructureManager  # noqa: E402

CANDIDATES = ROOT / "benchmarks" / "data" / "benchmark_candidates.csv"
STRATA = (("short", 0, 150), ("medium", 150, 300), ("long", 300, 10**9))

RESOLVED_COLUMNS = [
    "accession", "protein_name", "organism", "functional_class",
    "membrane_associated", "length", "stratum", "resolved_at",
]
RESULT_COLUMNS = [
    "accession", "protein_name", "organism", "functional_class", "stratum",
    "length", "provider", "status", "identifier", "mapping_level",
    "identity_percent", "coverage_percent", "accession_match", "confidence",
    "confidence_kind", "attempts", "recovered", "from_cache",
    "response_time_s", "error",
]


def stratum_for(length: int) -> str:
    for name, low, high in STRATA:
        if low <= length < high:
            return name
    return "long"  # pragma: no cover - unreachable given the bands above


def load_candidates(limit: int | None) -> list[dict]:
    with open(CANDIDATES, encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return rows[:limit] if limit else rows


def resolve(manager: StructureManager, candidates: list[dict], refresh: bool) -> tuple[list[dict], list[dict]]:
    """Fetch canonical sequences; return (resolved, unresolved)."""
    resolved, unresolved = [], []
    for index, row in enumerate(candidates, 1):
        accession = row["accession"].strip()
        print(f"  [{index}/{len(candidates)}] resolving {accession}", file=sys.stderr)
        try:
            entry = manager.uniprot.fetch_entry(accession, refresh=refresh)
        except BioSeqInsightError as exc:
            unresolved.append({**row, "reason": str(exc)})
            continue
        if not entry.sequence:
            unresolved.append({**row, "reason": "UniProt returned no sequence"})
            continue
        resolved.append(
            {
                **row,
                "sequence": entry.sequence,
                "length": entry.length,
                "stratum": stratum_for(entry.length),
                "uniprot_name": entry.name,
                "resolved_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            }
        )
    return resolved, unresolved


def benchmark(manager: StructureManager, resolved: list[dict], refresh: bool) -> list[dict]:
    rows: list[dict] = []
    for index, protein in enumerate(resolved, 1):
        accession = protein["accession"]
        print(
            f"  [{index}/{len(resolved)}] querying {accession} "
            f"({protein['length']} aa, {protein['stratum']})",
            file=sys.stderr,
        )
        outcome = manager.retrieve(accession, stop_on_exact=False, refresh=refresh)
        for result in outcome.results:
            validation = result.validation
            rows.append(
                {
                    "accession": accession,
                    "protein_name": protein["protein_name"],
                    "organism": protein["organism"],
                    "functional_class": protein["functional_class"],
                    "stratum": protein["stratum"],
                    "length": protein["length"],
                    "provider": result.source,
                    "status": result.status.value,
                    "identifier": result.identifier or "",
                    "mapping_level": result.mapping_level.value,
                    "identity_percent": _round(validation.identity_percent if validation else None),
                    "coverage_percent": _round(validation.coverage_percent if validation else None),
                    "accession_match": (validation.accession_match if validation else None),
                    "confidence": _round(result.confidence),
                    "confidence_kind": result.confidence_kind or "",
                    "attempts": len(result.attempts),
                    "recovered": result.recovered,
                    "from_cache": result.from_cache,
                    "response_time_s": _round(result.response_time_s, 3),
                    "error": (result.error or "").replace("\n", " ")[:200],
                }
            )
    return rows


def _round(value, digits: int = 2):
    return round(value, digits) if isinstance(value, (int, float)) else ""


def summarise(rows: list[dict], resolved: list[dict], unresolved: list[dict]) -> dict:
    providers = sorted({row["provider"] for row in rows})
    per_provider = {}
    for provider in providers:
        subset = [r for r in rows if r["provider"] == provider]
        ok = [r for r in subset if r["status"] == "ok"]
        times = [r["response_time_s"] for r in ok if isinstance(r["response_time_s"], (int, float))]
        not_found = sum(1 for r in subset if r["status"] == "not_found")
        timeouts = sum(1 for r in subset if r["status"] == "timeout")
        service_errors = sum(1 for r in subset if r["status"] == "service_error")
        invalid_input = sum(1 for r in subset if r["status"] == "invalid_input")
        skipped = sum(1 for r in subset if r["status"] == "skipped")
        # Every request lands in exactly one of these categories. Reporting
        # invalid_input and skipped separately -- rather than letting them
        # vanish from the breakdown -- matters here specifically: ESM Atlas
        # refuses sequences over its length limit before making any HTTP
        # request, so a run over this candidate set legitimately shows
        # dozens of requests with none of the other four statuses. Without
        # this field, "requests" and "retrieved" silently fail to reconcile
        # and a reader has no way to tell a length-limit rejection apart
        # from an unexplained gap in the accounting.
        accounted = len(ok) + not_found + timeouts + service_errors + invalid_input + skipped
        # requests actually sent over the wire, excluding client-side
        # refusals (invalid_input) and provider skips (offline mode). This
        # is the fairer "reliability" figure: success_rate_percent below
        # blends genuine service failures, genuine coverage gaps, AND
        # sequences that were never attempted (e.g. ESM Atlas's length
        # limit) into one denominator, which understates a provider that
        # simply cannot be asked about some candidates by design. Report
        # both rather than picking one, since they answer different
        # questions: "how often does this candidate set get a structure
        # from this provider" versus "how often does an attempted request
        # to this provider succeed".
        attempted = len(subset) - invalid_input - skipped
        per_provider[provider] = {
            "requests": len(subset),
            "retrieved": len(ok),
            "success_rate_percent": round(100.0 * len(ok) / len(subset), 1) if subset else 0.0,
            "attempted_requests": attempted,
            "attempted_success_rate_percent": (
                round(100.0 * len(ok) / attempted, 1) if attempted else None
            ),
            "not_found": not_found,
            "timeouts": timeouts,
            "service_errors": service_errors,
            "invalid_input": invalid_input,
            "skipped": skipped,
            "recovered_after_failure": sum(1 for r in subset if r["recovered"]),
            "median_response_time_s": round(_median(times), 3) if times else None,
            "exact_matches": sum(1 for r in subset if r["mapping_level"] == "M4"),
            "accounting_check": "ok" if accounted == len(subset) else (
                f"MISMATCH: {accounted} accounted for out of {len(subset)} requests"
            ),
        }

    # Best mapping achieved per protein, across all providers.
    # Every resolved candidate must appear in the distribution, including
    # ones that got nothing from any provider. Seeding every accession to
    # "M0" up front is what makes that true: the loop below only ever
    # *raises* a protein's level, so without this seed, a protein that
    # never gets past M0 anywhere is simply absent from `best` afterwards
    # rather than present with the value "M0" -- and an entry that was
    # never written is not counted by anything that only iterates
    # `best.values()`. A live run once produced a table whose levels summed
    # to one less than the resolved count for exactly this reason: one
    # protein had no usable structure from any resource and vanished from
    # the table instead of being tallied under M0.
    best: dict[str, str] = {p["accession"]: "M0" for p in resolved}
    for row in rows:
        current = best.get(row["accession"], "M0")
        if int(row["mapping_level"][1:]) > int(current[1:]):
            best[row["accession"]] = row["mapping_level"]
    distribution = {level.value: 0 for level in MappingLevel}
    for value in best.values():
        distribution[value] += 1
    assert sum(distribution.values()) == len(resolved), (
        f"best_mapping_distribution sums to {sum(distribution.values())}, "
        f"not {len(resolved)} resolved candidates -- this should be impossible "
        "now that every accession is seeded into `best`; please report this."
    )

    by_stratum: dict[str, dict] = {}
    for name, _, _ in STRATA:
        members = [p["accession"] for p in resolved if p["stratum"] == name]
        by_stratum[name] = {
            "proteins": len(members),
            "exact_matches": sum(1 for a in members if best.get(a) == "M4"),
        }

    return {
        "software_version": __version__,
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "candidates": len(resolved) + len(unresolved),
        "resolved": len(resolved),
        "unresolved": len(unresolved),
        "unresolved_accessions": [u["accession"] for u in unresolved],
        "strata": by_stratum,
        "per_provider": per_provider,
        "best_mapping_distribution": distribution,
        "exact_mapping_rate_percent": (
            round(100.0 * distribution["M4"] / len(resolved), 1) if resolved else 0.0
        ),
        "any_structure_rate_percent": (
            round(100.0 * sum(v for k, v in distribution.items() if k != "M0") / len(resolved), 1)
            if resolved
            else 0.0
        ),
    }


def _median(values: list[float]) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def write_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="benchmarks/results/structure_benchmark")
    parser.add_argument("--limit", type=int, help="only use the first N candidates")
    parser.add_argument("--refresh", action="store_true", help="ignore cached responses")
    parser.add_argument("--providers", help="comma-separated subset of rcsb,alphafold,esmatlas")
    parser.add_argument("--config", help="settings file to use")
    args = parser.parse_args(argv)

    settings = Settings.load(path=args.config)
    if args.providers:
        settings.fallback_order = tuple(p.strip() for p in args.providers.split(","))
    settings.validate()
    settings.ensure_directories()
    configure_logging(log_file=settings.log_file, level="INFO", console=False)

    if settings.offline:
        print("This benchmark requires network access; offline mode is enabled.", file=sys.stderr)
        return 2

    manager = StructureManager(settings)
    candidates = load_candidates(args.limit)
    print(f"Resolving {len(candidates)} candidate accessions against UniProt...", file=sys.stderr)
    started = time.perf_counter()
    resolved, unresolved = resolve(manager, candidates, args.refresh)

    if not resolved:
        print("No accession resolved; aborting.", file=sys.stderr)
        return 3
    print(
        f"Resolved {len(resolved)}/{len(candidates)}. Querying structural resources...",
        file=sys.stderr,
    )
    rows = benchmark(manager, resolved, args.refresh)
    elapsed = time.perf_counter() - started

    summary = summarise(rows, resolved, unresolved)
    summary["wall_clock_s"] = round(elapsed, 1)
    summary["settings"] = settings.to_dict()
    summary["cache"] = manager.cache_stats()

    out = Path(args.out)
    write_csv(out.with_name(out.name + "_resolved.csv"), resolved, RESOLVED_COLUMNS)
    write_csv(out.with_name(out.name + "_results.csv"), rows, RESULT_COLUMNS)
    out.with_name(out.name + "_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(json.dumps(
        {k: v for k, v in summary.items() if k not in ("settings", "cache")}, indent=2
    ))
    print(f"\nWritten to {out.parent}/", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
