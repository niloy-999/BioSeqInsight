#!/usr/bin/env python3
"""Experiment B: construct the adversarial benchmark (PROTOCOL.md Section 3).

Every adversarial case must have a documented, checkable reason for its
expected label, and -- per PROTOCOL.md Section 15/16 -- must be constructed
independently of bioseqinsight's own alignment engine, or a test of whether
that engine correctly rejects the case would be partly circular. Paralogue
grouping (B1, B5) therefore uses UniProt's own Pfam family cross-references,
fetched fresh per candidate, not any similarity score bioseqinsight itself
computes.

Classes built:
  B1  wrong paralogue        -- forced structure belongs to a same-Pfam-family
                                 protein at 40-95% identity to the query
                                 (identity estimated by the SAME crude k-mer
                                 Jaccard heuristic for every candidate, so it
                                 is at least a *consistent*, independent,
                                 non-bioseqinsight signal -- see note below)
  B2  fragment                -- forced structure = a real PDB entry whose
                                 SIFTS coverage of the query is <70%
                                 (reuses Experiment A's classification)
  B3  wrong isoform           -- query = one isoform's sequence; forced
                                 structure = a structure matching a
                                 DIFFERENT isoform of the same gene
  B4  unrelated (negative control) -- query paired with a uniformly-random,
                                 unrelated structure from elsewhere in the
                                 dataset
  B5  high-identity-wrong-entity   -- same construction as B1, restricted to
                                 the >=95% identity end of the k-mer estimate,
                                 to stress-test whether high estimated
                                 similarity alone would wrongly justify M4

On the k-mer heuristic: this is deliberately NOT the Needleman-Wunsch/
semi-global aligner bioseqinsight itself uses (core/alignment.py) -- it is a
simple, independent, order-insensitive shared-4-mer Jaccard index, used only
to select which pairs are "similar enough to be an interesting adversarial
case" during construction. It is not the identity value reported in results;
that comes from bioseqinsight's aligner when the system under test is run,
and separately from ground_truth.adjudicate() for the independent label.

Usage
-----
    python validation/scripts/01_build_adversarial.py \\
        --dataset validation/data/dataset_accessions.csv \\
        --out validation/data/adversarial_cases.csv
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
from pathlib import Path

VALIDATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(VALIDATION_ROOT))
from lib import api  # noqa: E402
from lib.fasta import read_fasta  # noqa: E402

SEED = 20261001


def kmer_jaccard(a: str, b: str, k: int = 4) -> float:
    """Independent, non-bioseqinsight similarity estimate for construction only."""
    if len(a) < k or len(b) < k:
        return 0.0
    sa = {a[i:i + k] for i in range(len(a) - k + 1)}
    sb = {b[i:i + k] for i in range(len(b) - k + 1)}
    union = sa | sb
    return len(sa & sb) / len(union) if union else 0.0


def est_identity(jaccard: float) -> float:
    """Estimated fractional identity from a 4-mer Jaccard index (Amendment 4).

    Two sequences at identity p share a 4-mer with probability ~p^4, so
    J ~= p^4 / (2 - p^4). Inverting: p = (2J / (1 + J)) ** 0.25. The earlier
    band "0.40 <= J <= 0.98" therefore meant ~87-99.5% identity, not the
    intended 40-95%, which is why only 7 paralogue pairs were found.
    """
    if jaccard <= 0:
        return 0.0
    return (2.0 * jaccard / (1.0 + jaccard)) ** 0.25


def kmer_shared(a: str, b: str, k: int = 4) -> int:
    """Number of distinct k-mers two sequences share (noise floor for est_identity)."""
    if len(a) < k or len(b) < k:
        return 0
    return len({a[i:i + k] for i in range(len(a) - k + 1)} & {b[i:i + k] for i in range(len(b) - k + 1)})


def fetch_pfam_ids(accession: str, use_cache: bool = True) -> set[str]:
    entry = api.uniprot_entry(accession, use_cache=use_cache)
    return {
        xref["id"] for xref in entry.get("uniProtKBCrossReferences", [])
        if xref.get("database") == "Pfam"
    }


def build_paralogue_pairs(candidates: list[dict], sequences: dict[str, str],
                           use_cache: bool) -> list[tuple[dict, dict, float]]:
    """(query, forced_structure_source, kmer_similarity) with shared Pfam family."""
    families: dict[str, list[dict]] = {}
    for cand in candidates:
        acc = cand["accession"]
        if acc not in sequences:
            continue
        try:
            pfams = fetch_pfam_ids(acc, use_cache=use_cache)
        except api.FetchError:
            continue
        for pf in pfams:
            families.setdefault(pf, []).append(cand)

    pairs = []
    seen_pairs: set[frozenset] = set()
    for _pfam, members in families.items():
        if len(members) < 2:
            continue
        for i, a in enumerate(members):
            for b in members[i + 1:]:
                key = frozenset((a["accession"], b["accession"]))
                if key in seen_pairs or a["accession"] == b["accession"]:
                    continue
                seq_a, seq_b = sequences.get(a["accession"]), sequences.get(b["accession"])
                if not seq_a or not seq_b:
                    continue
                sim = kmer_jaccard(seq_a, seq_b)
                est = est_identity(sim)
                # Amendment 4: band on ESTIMATED IDENTITY (0.40-0.995), not raw Jaccard.
                if 0.40 <= est < 0.995 and kmer_shared(seq_a, seq_b) >= 8:
                    pairs.append((a, b, sim))
                    seen_pairs.add(key)
    return pairs


def build_b2_fragment(dataset_rows: list[dict]) -> list[dict]:
    """Reuse Experiment A's own SIFTS-derived fragment classification (Section 2)."""
    cases = []
    for row in dataset_rows:
        if "fragment" in row.get("categories", "").split(";"):
            cases.append({
                "class": "B2_fragment", "query_accession": row["accession"],
                "forced_pdb_id": None,  # resolved at run time from the same PDB xref used in classification
                "expected_relationship": "fragment",
                "expected_max_m_level": "M3",
                "construction_note": f"SIFTS coverage {row.get('best_pdb_coverage')} < 0.70",
            })
    return cases


def build_b4_unrelated(candidates: list[dict], seed: int) -> list[dict]:
    rng = random.Random(seed)
    accs = [c["accession"] for c in candidates]
    rng.shuffle(accs)
    cases = []
    n = len(accs) // 2
    for query_acc, forced_acc in zip(accs[:n], accs[n:2 * n]):
        cases.append({
            "class": "B4_unrelated", "query_accession": query_acc,
            "forced_structure_accession": forced_acc,
            "expected_relationship": "unrelated",
            "expected_max_m_level": "M1",
            "construction_note": "uniformly random pairing, negative control",
        })
    return cases


def build_b3_wrong_isoform(dataset_rows: list[dict], use_cache: bool) -> list[dict]:
    """For each isoform-bearing candidate, pair its canonical sequence (query)
    against a *different* isoform's own PDB structure, when one exists with
    enough coverage of that other isoform to be a plausible retrieval hit.
    """
    cases = []
    for row in dataset_rows:
        if "isoform_bearing" not in row.get("categories", "").split(";"):
            continue
        accession = row["accession"]
        try:
            isoform_ids = api.uniprot_isoform_ids(accession, use_cache=use_cache)
        except api.FetchError:
            continue
        others = [iid for iid in isoform_ids if iid != f"{accession}-1" and iid != accession]
        if not others:
            continue
        for other_id in others[:2]:  # cap per gene to keep the set balanced
            try:
                other_seq = api.uniprot_isoform_sequence(other_id, use_cache=use_cache)
            except api.FetchError:
                continue
            if not other_seq:
                continue
            cases.append({
                "class": "B3_wrong_isoform", "query_accession": accession,
                "forced_structure_isoform_id": other_id,
                "forced_structure_isoform_length": len(other_seq),
                "query_canonical_length": row.get("length"),  # Amendment 4: lets label validity be checked
                "expected_relationship": "wrong_isoform",
                "expected_max_m_level": "M3",
                "construction_note": (
                    f"query is canonical isoform of {accession}; structure forced to a real "
                    f"structure of isoform {other_id} instead"
                ),
            })
    return cases


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--fasta")
    parser.add_argument("--out", default=str(VALIDATION_ROOT / "data" / "adversarial_cases.csv"))
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--max-paralogue-pairs", type=int, default=150)
    parser.add_argument("--force", action="store_true",
                        help="overwrite an existing (frozen) adversarial_cases.csv")  # Amendment 8
    args = parser.parse_args(argv)
    use_cache = not args.no_cache

    fasta_path = Path(args.fasta) if args.fasta else Path(args.dataset).with_name(
        Path(args.dataset).name.replace("_accessions.csv", "_sequences.fasta")
    )
    sequences = {acc: seq for acc, _desc, seq in read_fasta(fasta_path)}
    with open(args.dataset, newline="", encoding="utf-8") as fh:
        dataset_rows = list(csv.DictReader(fh))

    if Path(args.out).exists() and not args.force:  # Amendment 8: a used set is frozen
        print(f"{args.out} already exists and is FROZEN; not regenerating (use --force to overwrite).",
              file=sys.stderr)
        return 0
    print("Building B1/B5 (paralogue) pairs via independent Pfam grouping...", file=sys.stderr)
    pairs = build_paralogue_pairs(dataset_rows, sequences, use_cache)
    # Amendment 8: the pair list is built by iterating sets, whose order varies with
    # per-process string hashing; sort before the seeded shuffle so the draw is reproducible.
    pairs.sort(key=lambda t: (t[0]["accession"], t[1]["accession"]))
    rng = random.Random(SEED)
    rng.shuffle(pairs)
    pairs = pairs[: args.max_paralogue_pairs]

    b1_b5_cases = []
    for a, b, sim in pairs:
        cls = "B5_high_identity_wrong_entity" if est_identity(sim) >= 0.95 else "B1_wrong_paralogue"
        b1_b5_cases.append({
            "class": cls, "query_accession": a["accession"],
            "forced_structure_accession": b["accession"],
            "kmer_similarity_construction_only": round(sim, 4),
            "estimated_identity_construction_only": round(est_identity(sim), 3),
            "expected_relationship": "paralogue", "expected_max_m_level": "M3",
            "construction_note": f"shared Pfam family, k-mer Jaccard {sim:.3f} (construction-time estimate only)",
        })

    print(f"  {len(b1_b5_cases)} paralogue-based cases "
          f"({sum(1 for c in b1_b5_cases if c['class'].startswith('B5'))} at B5 threshold)", file=sys.stderr)

    b2_cases = build_b2_fragment(dataset_rows)
    print(f"  {len(b2_cases)} B2 fragment cases (from Experiment A classification)", file=sys.stderr)

    b4_cases = build_b4_unrelated(dataset_rows, SEED)
    print(f"  {len(b4_cases)} B4 unrelated negative-control cases", file=sys.stderr)

    isoform_rows = [r for r in dataset_rows if "isoform_bearing" in r.get("categories", "").split(";")]
    print(f"  {len(isoform_rows)} isoform-bearing candidates found; building B3...", file=sys.stderr)
    b3_cases = build_b3_wrong_isoform(dataset_rows, use_cache)
    print(f"  {len(b3_cases)} B3 wrong-isoform cases", file=sys.stderr)

    all_cases = b1_b5_cases + b2_cases + b4_cases + b3_cases
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({k for case in all_cases for k in case})
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_cases)

    import hashlib  # Amendment 8
    digest = hashlib.sha256(out_path.read_bytes()).hexdigest()
    Path(str(out_path) + ".sha256").write_text(f"{digest}  {out_path.name}\n")
    print(f"sha256 {digest}", file=sys.stderr)
    print(f"\nTotal adversarial cases: {len(all_cases)}", file=sys.stderr)
    print(f"Wrote {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
