#!/usr/bin/env python3
"""Experiment A: build the stratified 1,000-protein benchmark automatically.

This performs PROTOCOL.md Section 2 end to end with no manual curation:
queries UniProtKB directly for each length stratum, fetches full entries,
classifies each candidate into the six categories by calling RCSB and
PDBe SIFTS independently, draws a seeded stratified sample, and writes
both the frozen accession list and a ready-to-use FASTA file.

Why UniProt rather than NCBI/RefSeq as the source of accessions: every
resource `bioseqinsight` queries (RCSB's cross-reference check, AlphaFold
DB, the UniProt canonical-sequence fetch itself) is keyed by UniProt
accession, not a RefSeq or GenBank identifier. Building the benchmark from
NCBI accessions would require an extra, error-prone identifier-mapping step
before the software under test could even be invoked. InterPro is used
below only as an optional refinement for paralogue-family detection
(Experiment B), not as a primary sequence source.

Usage
-----
    python validation/scripts/00_build_dataset.py --target 100   # pilot
    python validation/scripts/00_build_dataset.py --target 1000  # main run

Everything is cached on disk (validation/lib/api.py), so re-running after an
interruption resumes quickly rather than re-downloading.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(VALIDATION_ROOT))
from lib import api  # noqa: E402
from lib.fasta import write_fasta  # noqa: E402

SEED = 20261001
STRATA = [
    ("short", 50, 150),
    ("medium", 150, 300),
    ("long", 300, 600),
    ("very_long", 600, 100000),
]
CATEGORIES = [
    "exact_full_length", "fragment", "partial_missing_residue",
    "isoform_bearing", "paralogue_family", "no_structure",
]
TARGET_PER_STRATUM_FRACTION = 0.25  # four strata, ~equal split


@dataclass
class Candidate:
    accession: str
    length: int
    stratum: str
    organism: str = ""
    protein_name: str = ""
    sequence: str = ""
    pdb_ids: list[str] = field(default_factory=list)
    has_isoforms: bool = False
    categories: list[str] = field(default_factory=list)
    pfam_ids: list[str] = field(default_factory=list)  # Amendment 2
    best_pdb_coverage: float | None = None   # fraction of canonical seq covered by best SEQRES
    has_unobserved_residues: bool | None = None
    resolved: bool = True
    note: str = ""

    def to_row(self) -> dict:
        return {
            "accession": self.accession, "length": self.length, "stratum": self.stratum,
            "organism": self.organism, "protein_name": self.protein_name,
            "pdb_ids": ";".join(self.pdb_ids), "has_isoforms": self.has_isoforms,
            "categories": ";".join(self.categories),
            "pfam_ids": ";".join(self.pfam_ids),  # Amendment 2
            "best_pdb_coverage": self.best_pdb_coverage,
            "has_unobserved_residues": self.has_unobserved_residues,
            "resolved": self.resolved, "note": self.note,
        }


def uniprot_query_for_stratum(low: int, high: int) -> str:
    """Documented UniProtKB REST query syntax (reviewed entries only, length range).

    If UniProt's query grammar has changed since this was written, this is
    the single line to fix; everything downstream classifies candidates from
    the fetched JSON rather than depending on further query-string tricks.
    """
    hi = high if high < 100000 else 100000
    return f"reviewed:true AND length:[{low} TO {hi}]"


def fetch_stratum_pool(stratum: str, low: int, high: int, pool_size: int,
                        use_cache: bool = True) -> list[dict]:
    """Page through UniProt search results for one stratum up to pool_size entries."""
    query = uniprot_query_for_stratum(low, high)
    print(f"  [UniProt] querying stratum '{stratum}': {query!r}", file=sys.stderr)
    results: list[dict] = []
    cursor = None
    release = None
    while len(results) < pool_size:
        batch, cursor, release = api.uniprot_search(query, size=500, cursor=cursor, use_cache=use_cache)
        if not batch:
            break
        results.extend(batch)
        print(f"    fetched {len(results)} so far (UniProt release {release})", file=sys.stderr)
        if not cursor:
            break
    return results[:pool_size]


SAMPLING_NOTES: list[str] = []  # Amendment 3: surfaced in the output metadata


def sampling_diagnostics(chosen: list) -> dict:
    """How concentrated is the final sample? Reported so bias is visible."""
    import collections
    organisms = collections.Counter(c.organism for c in chosen)
    prefixes = collections.Counter(c.accession[:2] for c in chosen)
    top_org, top_n = organisms.most_common(1)[0] if organisms else ("", 0)
    return {
        "sampling_notes": list(SAMPLING_NOTES),
        "distinct_organisms": len(organisms),
        "top_organism": top_org,
        "top_organism_fraction": round(top_n / len(chosen), 3) if chosen else None,
        "top_accession_prefixes": prefixes.most_common(5),
    }


def fetch_stratum_pool_random(stratum: str, low: int, high: int, pool_size: int,
                               use_cache: bool = True) -> list[dict]:
    """Amendment 3: a seeded RANDOM pool for one stratum.

    Streams every accession in the length range, samples ``pool_size`` of
    them with a per-stratum seeded RNG, then fetches those entries. The
    earlier fetch_stratum_pool() took the first N search results, which the
    pilot showed were 100% Homo sapiens and alphabetically clustered.
    """
    query = uniprot_query_for_stratum(low, high)
    try:
        accessions, release = api.uniprot_stream_accessions(query, use_cache=use_cache)
    except api.FetchError as exc:
        note = f"{stratum}: UniProt stream failed ({exc}); fell back to first-N paging, SAMPLE IS BIASED"
        SAMPLING_NOTES.append(note)
        print(f"  ! {note}", file=sys.stderr)
        return fetch_stratum_pool(stratum, low, high, pool_size, use_cache=use_cache)

    universe = sorted(set(accessions))
    rng = random.Random(f"{SEED}-{stratum}")
    picked = rng.sample(universe, min(pool_size, len(universe)))
    SAMPLING_NOTES.append(
        f"{stratum}: random sample of {len(picked)} from {len(universe)} reviewed entries "
        f"(UniProt release {release})"
    )
    print(f"  [UniProt] stratum '{stratum}': {len(picked)} random of {len(universe)} entries", file=sys.stderr)

    entries: list[dict] = []
    for i, accession in enumerate(picked, 1):
        try:
            entries.append(api.uniprot_entry(accession, use_cache=use_cache))
        except api.FetchError as exc:
            print(f"    ! skipping {accession}: {exc}", file=sys.stderr)
        if i % 100 == 0:
            print(f"    fetched {i}/{len(picked)} entries", file=sys.stderr)
    return entries


def classify_candidate(entry_summary: dict, use_cache: bool = True) -> Candidate:
    """Classify one UniProt search-result entry into A's categories (Section 2).

    Calls RCSB and SIFTS independently -- this is dataset construction, not
    the system under test, and PROTOCOL.md Section 5 requires ground truth
    and categorisation to be independent of bioseqinsight's own alignment.
    """
    accession = entry_summary["primaryAccession"]
    length = entry_summary["sequence"]["length"]
    organism = entry_summary.get("organism", {}).get("scientificName", "")
    protein_name = (
        entry_summary.get("proteinDescription", {})
        .get("recommendedName", {}).get("fullName", {}).get("value", "")
    )
    pdb_ids = [
        xref["id"] for xref in entry_summary.get("uniProtKBCrossReferences", [])
        if xref.get("database") == "PDB"
    ]
    pfam_ids = [  # Amendment 2
        xref["id"] for xref in entry_summary.get("uniProtKBCrossReferences", [])
        if xref.get("database") == "Pfam"
    ]
    has_isoforms = bool(entry_summary.get("comments")) and any(
        c.get("commentType") == "ALTERNATIVE PRODUCTS" for c in entry_summary.get("comments", [])
    )

    cand = Candidate(
        accession=accession, length=length, stratum="",  # filled by caller
        organism=organism, protein_name=protein_name,
        pdb_ids=pdb_ids, has_isoforms=has_isoforms, pfam_ids=pfam_ids,  # Amendment 2
        sequence=entry_summary.get("sequence", {}).get("value", ""),
    )

    if has_isoforms:
        cand.categories.append("isoform_bearing")

    if not pdb_ids:
        cand.categories.append("no_structure")
        return cand

    # Independent ground truth for coverage/completeness: SIFTS residue mapping.
    try:
        sifts = api.sifts_uniprot_mappings(accession, use_cache=use_cache)
    except api.FetchError as exc:
        cand.note = f"SIFTS fetch failed: {exc}"
        sifts = {}

    best_coverage = 0.0
    # best_structures returns a FLAT LIST directly under the accession key
    # (PROTOCOL.md Amendment 1: confirmed against two independent client
    # implementations). No "PDB" sub-key, no per-entity "mappings" sub-list:
    # each list entry IS one PDB structure's mapping, carrying
    # pdb_id/chain_id/coverage/unp_start/unp_end/resolution directly. Try
    # the accession as given, then uppercased, in case PDBe's key casing
    # ever differs from UniProt's.
    entries = []
    if sifts:
        entries = sifts.get(accession) or sifts.get(accession.upper()) or []
    for entry in entries:
        if "coverage" in entry and entry["coverage"] is not None:
            # PDBe computes this itself; prefer it over a locally-derived value.
            coverage = float(entry["coverage"])
            if coverage > 1.0:  # defensively handle a 0-100 scale if PDBe ever returns one
                coverage /= 100.0
        else:
            start, end = entry.get("unp_start"), entry.get("unp_end")
            if start is None or end is None or not length:
                continue
            coverage = (end - start + 1) / length
        best_coverage = max(best_coverage, coverage)
    # Per-residue unresolved/gap information is not part of this endpoint's
    # response (see the correction note in lib/api.py); left unpopulated
    # rather than guessed from a field that was never confirmed to exist.
    any_unobserved = None

    cand.best_pdb_coverage = round(best_coverage, 4) if entries else None
    cand.has_unobserved_residues = any_unobserved

    if best_coverage >= 0.90:
        cand.categories.append("exact_full_length")
    elif entries and best_coverage < 0.70:
        cand.categories.append("fragment")
    if not entries:
        # PDB cross-reference existed in UniProt but SIFTS has no mapping
        # (rare; usually a very recent deposition). Record and keep for
        # manual inspection rather than silently guessing a category.
        cand.note = (cand.note + "; " if cand.note else "") + "PDB xref present but no SIFTS mapping found"

    return cand


def stratified_sample(pool: dict[str, list[Candidate]], target: int, seed: int) -> list[Candidate]:
    """Seeded stratified draw: ~target/4 per length stratum, spread across categories."""
    rng = random.Random(seed)
    per_stratum = target // len(STRATA)
    chosen: list[Candidate] = []
    for stratum, _low, _high in STRATA:
        candidates = pool.get(stratum, [])
        by_category: dict[str, list[Candidate]] = {c: [] for c in CATEGORIES}
        for cand in candidates:
            cats = cand.categories or ["no_structure"]
            for cat in cats:
                by_category[cat].append(cand)
        picked_ids: set[str] = set()
        per_category = max(1, per_stratum // len(CATEGORIES))
        for cat in CATEGORIES:
            bucket = [c for c in by_category[cat] if c.accession not in picked_ids]
            rng.shuffle(bucket)
            for cand in bucket[:per_category]:
                picked_ids.add(cand.accession)
                chosen.append(cand)
        # top up the stratum to its quota from whatever is left, if categories
        # were sparse (expected for e.g. "isoform_bearing" in the short stratum)
        remaining = [c for c in candidates if c.accession not in picked_ids]
        rng.shuffle(remaining)
        have = sum(1 for c in chosen if c.stratum == stratum)
        need = per_stratum - have
        for cand in remaining[: max(0, need)]:
            picked_ids.add(cand.accession)
            chosen.append(cand)
    return chosen


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--target", type=int, default=1000, help="total dataset size")
    parser.add_argument("--pool-multiplier", type=float, default=3.0,
                         help="fetch this many candidates per stratum before classifying/sampling")
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--out-prefix", default=str(VALIDATION_ROOT / "data" / "dataset"))
    args = parser.parse_args(argv)
    use_cache = not args.no_cache

    started = time.time()
    pool_size_per_stratum = int(args.target / len(STRATA) * args.pool_multiplier)
    pool: dict[str, list[Candidate]] = {s: [] for s, _, _ in STRATA}
    release_seen = None

    for stratum, low, high in STRATA:
        raw = fetch_stratum_pool_random(stratum, low, high, pool_size_per_stratum, use_cache=use_cache)
        print(f"  classifying {len(raw)} candidates in stratum '{stratum}'...", file=sys.stderr)
        for i, entry in enumerate(raw, 1):
            try:
                cand = classify_candidate(entry, use_cache=use_cache)
            except Exception as exc:  # noqa: BLE001 - one bad record must not abort the run
                print(f"    ! skipping {entry.get('primaryAccession', '?')}: {exc}", file=sys.stderr)
                continue
            cand.stratum = stratum
            pool[stratum].append(cand)
            if i % 50 == 0:
                print(f"    {i}/{len(raw)} classified", file=sys.stderr)

    # Amendment 2: tag paralogue_family across the WHOLE pool -- a
    # cross-candidate, group-level property (a protein only counts as a
    # paralogue-family member if >=1 OTHER candidate anywhere in the pool
    # shares a Pfam ID with it), using Pfam IDs already fetched as part of
    # the UniProt search result, with no extra network calls.
    pfam_groups: dict[str, list] = {}
    for stratum_candidates in pool.values():
        for cand in stratum_candidates:
            for pfam_id in cand.pfam_ids:
                pfam_groups.setdefault(pfam_id, []).append(cand)
    tagged = set()
    for pfam_id, members in pfam_groups.items():
        if len(members) >= 2:
            for cand in members:
                if "paralogue_family" not in cand.categories:
                    cand.categories.append("paralogue_family")
                tagged.add(cand.accession)
    print(f"  tagged {len(tagged)} candidates as paralogue_family "
          f"(sharing a Pfam ID with >=1 other candidate in the pool)", file=sys.stderr)

    chosen = stratified_sample(pool, args.target, SEED)
    print(f"\nSelected {len(chosen)} candidates (target {args.target}).", file=sys.stderr)

    accessions_path = Path(f"{args.out_prefix}_accessions.csv")
    with open(accessions_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(chosen[0].to_row().keys()) if chosen else [])
        writer.writeheader()
        for cand in chosen:
            writer.writerow(cand.to_row())

    fasta_path = Path(f"{args.out_prefix}_sequences.fasta")
    write_fasta(fasta_path, [(c.accession, f"{c.accession} {c.protein_name}", c.sequence) for c in chosen if c.sequence])

    category_counts: dict[str, int] = {c: 0 for c in CATEGORIES}
    for cand in chosen:
        for cat in (cand.categories or ["no_structure"]):
            category_counts[cat] += 1
    stratum_counts = {s: sum(1 for c in chosen if c.stratum == s) for s, _, _ in STRATA}

    metadata = {
        "seed": SEED,
        "target": args.target,
        "achieved": len(chosen),
        "elapsed_s": round(time.time() - started, 1),
        "stratum_counts": stratum_counts,
        "category_counts": category_counts,
        "sampling_diagnostics": sampling_diagnostics(chosen),
    }
    Path(f"{args.out_prefix}_metadata.json").write_text(json.dumps(metadata, indent=2))

    print(json.dumps(metadata, indent=2))
    print(f"\nWrote {accessions_path}\nWrote {fasta_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
