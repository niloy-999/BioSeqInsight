"""Thin HTTP clients for the data sources in PROTOCOL.md Section 1.

Deliberately independent of the installed ``bioseqinsight`` package's own
HTTP stack: this code is the thing that BUILDS the ground truth and the
candidate list, and Section 5 of the protocol requires ground truth to be
established independently of the system under test. Sharing a code path
with ``bioseqinsight`` here would not be circular in a statistical sense,
but keeping the two separate removes any doubt and makes this package
usable even against a ``bioseqinsight`` checkout that does not import
cleanly.

Everything here is stdlib ``urllib`` plus a simple on-disk JSON cache and a
token-bucket rate limiter, so nothing extra needs to be installed to run
the benchmark scripts.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

USER_AGENT = "BioSeqInsight-Validation/1.0 (research benchmark; contact via GitHub repo)"
CACHE_DIR = Path(os.environ.get("BSI_VALIDATION_CACHE", Path.home() / ".bioseqinsight" / "validation-cache"))


@dataclass
class RateLimiter:
    """Simple fixed-interval limiter: at most one request per ``interval_s``."""

    interval_s: float
    _last: float = field(default=0.0, init=False)

    def wait(self) -> None:
        elapsed = time.monotonic() - self._last
        if elapsed < self.interval_s:
            time.sleep(self.interval_s - elapsed)
        self._last = time.monotonic()


class FetchError(RuntimeError):
    """Raised after the retry budget for a single request is exhausted."""

    def __init__(self, url: str, status: int | None, detail: str):
        self.url = url
        self.status = status
        self.detail = detail
        super().__init__(f"{status or 'network error'} fetching {url}: {detail}")


def _cache_key(url: str) -> str:
    return hashlib.sha256(url.encode("utf-8")).hexdigest()


def _cache_path(namespace: str, url: str) -> Path:
    d = CACHE_DIR / namespace
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{_cache_key(url)}.json"


def cached_get_json(
    url: str,
    *,
    namespace: str,
    limiter: RateLimiter | None = None,
    max_attempts: int = 4,
    backoff_base_s: float = 1.0,
    timeout_s: float = 20.0,
    use_cache: bool = True,
) -> tuple[Any, dict]:
    """GET ``url``, parse JSON, cache on disk. Returns ``(payload, meta)``.

    ``meta`` carries ``{"from_cache": bool, "headers": {...}, "status": int}``
    so callers (in particular the UniProt release-header logic in
    PROTOCOL.md Section 1) can inspect response headers even on a cache hit,
    since headers are cached alongside the body.
    """
    path = _cache_path(namespace, url)
    if use_cache and path.is_file():
        try:
            cached = json.loads(path.read_text(encoding="utf-8"))
            return cached["body"], {**cached["meta"], "from_cache": True}
        except (json.JSONDecodeError, KeyError):
            pass  # corrupt cache entry: refetch

    last_error: Exception | None = None
    for attempt in range(1, max_attempts + 1):
        if limiter is not None:
            limiter.wait()
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
        try:
            with urllib.request.urlopen(request, timeout=timeout_s) as response:
                raw = response.read()
                headers = dict(response.headers.items())
                status = response.status
            body = json.loads(raw.decode("utf-8")) if raw else None
            meta = {"headers": headers, "status": status, "from_cache": False, "url": url}
            if use_cache:
                path.write_text(json.dumps({"body": body, "meta": meta}), encoding="utf-8")
            return body, meta
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                raise FetchError(url, 404, "not found") from exc
            if exc.code in (400, 401, 403):
                raise FetchError(url, exc.code, "client error, not retried") from exc
            last_error = exc
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last_error = exc
        if attempt < max_attempts:
            time.sleep(backoff_base_s * (2 ** (attempt - 1)))
    raise FetchError(url, getattr(last_error, "code", None), str(last_error))


# ---------------------------------------------------------------------------
# UniProt
# ---------------------------------------------------------------------------

_uniprot_limiter = RateLimiter(interval_s=0.34)  # ~3 req/s, polite default


def uniprot_search(query: str, *, size: int = 500, fields: str | None = None,
                    cursor: str | None = None, use_cache: bool = True) -> tuple[list[dict], str | None, str | None]:
    """Search UniProtKB. Returns ``(results, next_cursor, uniprot_release)``."""
    fields = fields or ("accession,id,protein_name,organism_name,length,"
                        "cc_alternative_products,xref_pdb,xref_pfam")  # Amendment 2: + xref_pfam
    url = (
        "https://rest.uniprot.org/uniprotkb/search"
        f"?query={urllib.parse.quote(query)}&format=json&size={size}&fields={fields}"
    )
    if cursor:
        url += f"&cursor={cursor}"
    body, meta = cached_get_json(url, namespace="uniprot_search", limiter=_uniprot_limiter, use_cache=use_cache)
    results = body.get("results", []) if body else []
    next_cursor = None
    link = meta["headers"].get("Link", "")
    if 'rel="next"' in link:
        for part in link.split(","):
            if 'rel="next"' in part and "cursor=" in part:
                next_cursor = part.split("cursor=")[1].split("&")[0].split(">")[0]
    release = meta["headers"].get("X-UniProt-Release")
    return results, next_cursor, release


def uniprot_entry(accession: str, *, use_cache: bool = True) -> dict:
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    body, _meta = cached_get_json(url, namespace="uniprot_entry", limiter=_uniprot_limiter, use_cache=use_cache)
    return body


# ---------------------------------------------------------------------------
# RCSB PDB
# ---------------------------------------------------------------------------

_rcsb_limiter = RateLimiter(interval_s=0.25)


def rcsb_polymer_entity(pdb_id: str, entity_id: str = "1", *, use_cache: bool = True) -> dict:
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
    body, _meta = cached_get_json(url, namespace="rcsb_entity", limiter=_rcsb_limiter, use_cache=use_cache)
    return body


def rcsb_entry(pdb_id: str, *, use_cache: bool = True) -> dict:
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    body, _meta = cached_get_json(url, namespace="rcsb_entry", limiter=_rcsb_limiter, use_cache=use_cache)
    return body


# ---------------------------------------------------------------------------
# PDBe SIFTS (independent cross-reference, used as ground truth per Section 5)
# ---------------------------------------------------------------------------

_sifts_limiter = RateLimiter(interval_s=0.34)


def sifts_uniprot_mappings(accession: str, *, use_cache: bool = True) -> dict:
    """PDB structures mapped to this UniProt accession, via PDBe's `best_structures`.

    CORRECTED (validation/PROTOCOL.md Amendment 1): an earlier version of
    this function called ``/mappings/uniprot/{accession}``, which is the
    *reverse* lookup -- documented usage of that path takes a **PDB ID** and
    returns the UniProt accession(s) it maps to (e.g.
    ``/mappings/uniprot/1ivv``), not the other way round. Called with a
    UniProt accession, it does not return the per-structure coverage data
    this function's callers need, which is exactly why the pilot run's
    exact_full_length/fragment/partial_missing_residue categories came back
    empty regardless of whether a candidate actually had a PDB structure.

    The correct UniProt-accession-keyed endpoint is PDBe's own
    ``best_structures`` API: ``/mappings/best_structures/{accession}``,
    confirmed against two independent published client implementations in
    addition to this project's own documentation trail. Its response shape
    is ``{accession: [entry, entry, ...]}`` -- a **flat list** directly
    under the accession key, each entry carrying (at least) ``pdb_id``,
    ``chain_id``, ``coverage`` (already computed by PDBe, 0-1 scale),
    ``resolution``, ``unp_start`` and ``unp_end``. There is no "PDB"
    sub-key to nest under, which was this function's second bug.

    Per-residue gap (unresolved-residue) information is NOT part of this
    response and is not fabricated here; ``partial_missing_residue``
    classification in ``00_build_dataset.py`` is consequently left
    unpopulated in this version rather than guessed from an unverified
    field -- see validation/README.md Section 6.
    """
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{accession}"
    try:
        body, _meta = cached_get_json(url, namespace="sifts_best_structures",
                                       limiter=_sifts_limiter, use_cache=use_cache)
    except FetchError as exc:
        if exc.status == 404:
            return {}
        raise
    return body or {}


def rcsb_sequence_search(sequence: str, *, identity_cutoff: float | None = None,
                          evalue_cutoff: float = 1.0, max_hits: int = 10,
                          use_cache: bool = True) -> list[dict]:
    """RCSB's own sequence-search API (POST, so not routed through cached_get_json's
    GET-only cache; a lightweight local cache keyed on the request body is used
    instead so repeated runs against the same sequence do not re-query RCSB).

    ``identity_cutoff`` is RCSB's own request-time filter (``identity_cutoff``
    search-node parameter, 0-1). Passing ``None`` performs an unfiltered
    top-hit search (Method 1 in PROTOCOL.md Section 4); passing a cutoff
    (e.g. 0.9, bioseqinsight's own default ``rcsb_search_identity_cutoff``)
    performs RCSB's own filtered search (Method 2). Both are RCSB's own
    ranking; the difference is only whether RCSB's own filter was applied
    at request time, not any post-hoc processing by this code or by
    bioseqinsight.

    Returns the raw list of hits (``[{"identifier": "1ABC_1", "score": ...}, ...]``),
    highest score first, as RCSB returns them.
    """
    node: dict = {
        "type": "terminal",
        "service": "sequence",
        "parameters": {
            "evalue_cutoff": evalue_cutoff,
            "sequence_type": "protein",
            "value": sequence,
        },
    }
    if identity_cutoff is not None:
        node["parameters"]["identity_cutoff"] = identity_cutoff
    payload = {
        "query": node,
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": 0, "rows": max_hits}},
    }
    body_key = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    cache_file = _cache_path("rcsb_search", body_key)
    if use_cache and cache_file.is_file():
        try:
            return json.loads(cache_file.read_text())["result_set"]
        except (json.JSONDecodeError, KeyError):
            pass

    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    data = json.dumps(payload).encode("utf-8")
    last_error: Exception | None = None
    for attempt in range(1, 4):
        _rcsb_limiter.wait()
        request = urllib.request.Request(
            url, data=data, method="POST",
            headers={"User-Agent": USER_AGENT, "Content-Type": "application/json"},
        )
        try:
            with urllib.request.urlopen(request, timeout=25.0) as response:
                raw = response.read().decode("utf-8")
            # Amendment 4: RCSB answers a search with no hits as HTTP 204 and an
            # EMPTY body, which urllib returns normally (it is not an HTTPError).
            if not raw.strip():
                result_set = []
            else:
                try:
                    body = json.loads(raw)
                except ValueError as decode_error:
                    # URLError is retried by the handler below, then raised as FetchError.
                    raise urllib.error.URLError(f"invalid JSON from RCSB: {decode_error}")
                result_set = body.get("result_set", [])
            if use_cache:
                cache_file.write_text(json.dumps({"result_set": result_set}))
            return result_set
        except urllib.error.HTTPError as exc:
            if exc.code == 204 or exc.code == 404:
                return []  # RCSB returns 204/404-shaped "no hits" for some empty searches
            last_error = exc
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last_error = exc
        time.sleep(1.5 * attempt)
    raise FetchError(url, getattr(last_error, "code", None), str(last_error))


def fetch_fasta_text(url: str, *, namespace: str, limiter: RateLimiter | None = None,
                      max_attempts: int = 3, use_cache: bool = True) -> str:
    """GET a plain-text (non-JSON) response, e.g. a UniProt isoform FASTA record.

    UniProt serves individual isoforms (accession suffix ``-2``, ``-3``, ...)
    only via the plain FASTA/text endpoints, not the JSON entry endpoint used
    by :func:`uniprot_entry`, hence a separate text-mode fetch path.
    """
    path = _cache_path(namespace, url)
    if use_cache and path.is_file():
        cached = json.loads(path.read_text(encoding="utf-8"))
        return cached["text"]
    last_error: Exception | None = None
    for attempt in range(1, max_attempts + 1):
        if limiter is not None:
            limiter.wait()
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=20.0) as response:
                text = response.read().decode("utf-8")
            if use_cache:
                path.write_text(json.dumps({"text": text}), encoding="utf-8")
            return text
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                raise FetchError(url, 404, "not found") from exc
            last_error = exc
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last_error = exc
        time.sleep(1.0 * attempt)
    raise FetchError(url, getattr(last_error, "code", None), str(last_error))


def uniprot_isoform_sequence(isoform_id: str, *, use_cache: bool = True) -> str:
    """Fetch one isoform's sequence (e.g. accession ``P04637-2``) as plain sequence."""
    url = f"https://rest.uniprot.org/uniprotkb/{isoform_id}.fasta"
    text = fetch_fasta_text(url, namespace="uniprot_isoform", limiter=_uniprot_limiter, use_cache=use_cache)
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith(">")]
    return "".join(lines)


def uniprot_isoform_ids(accession: str, *, use_cache: bool = True) -> list[str]:
    """Isoform accessions (e.g. ``["P04637-1", "P04637-2"]``) from the full entry."""
    entry = uniprot_entry(accession, use_cache=use_cache)
    ids: list[str] = []
    for comment in entry.get("comments", []):
        if comment.get("commentType") != "ALTERNATIVE PRODUCTS":
            continue
        for iso in comment.get("isoforms", []):
            for iso_id in iso.get("isoformIds", []):
                ids.append(iso_id)
    return ids


def fetch_rcsb_pdb_file(pdb_id: str, *, use_cache: bool = True) -> str:
    """Raw coordinate file text, independent of bioseqinsight's own download path."""
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    return fetch_fasta_text(url, namespace="rcsb_pdb_file", limiter=_rcsb_limiter, use_cache=use_cache)


def uniprot_stream_accessions(query: str, *, use_cache: bool = True) -> tuple[list[str], str | None]:
    """Every accession matching ``query``, via UniProt's stream endpoint.

    Amendment 3. Verified live: ``format=tsv`` returns a header line
    ("Entry") followed by one accession per line, sorted alphabetically.
    Returns ``(accessions, uniprot_release)``.
    """
    url = ("https://rest.uniprot.org/uniprotkb/stream?query=" + urllib.parse.quote(query)
           + "&fields=accession&format=tsv")
    path = _cache_path("uniprot_stream", url)
    if use_cache and path.is_file():
        try:
            cached = json.loads(path.read_text(encoding="utf-8"))
            return cached["accessions"], cached.get("release")
        except (json.JSONDecodeError, KeyError):
            pass
    last_error: Exception | None = None
    for attempt in range(1, 4):
        _uniprot_limiter.wait()
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=180.0) as response:
                text = response.read().decode("utf-8")
                release = response.headers.get("X-UniProt-Release")
            lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
            accessions = [ln for ln in lines if ln != "Entry"]
            if not accessions:
                raise FetchError(url, None, "stream returned no accessions")
            if use_cache:
                path.write_text(json.dumps({"accessions": accessions, "release": release}), encoding="utf-8")
            return accessions, release
        except FetchError:
            raise
        except urllib.error.HTTPError as exc:
            last_error = exc
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last_error = exc
        time.sleep(2.0 * attempt)
    raise FetchError(url, getattr(last_error, "code", None), str(last_error))


__all__ = [
    "CACHE_DIR", "FetchError", "RateLimiter",
    "cached_get_json", "fetch_fasta_text", "fetch_rcsb_pdb_file", "uniprot_search", "uniprot_entry",
    "uniprot_isoform_sequence", "uniprot_isoform_ids",
    "rcsb_polymer_entity", "rcsb_entry", "rcsb_sequence_search", "sifts_uniprot_mappings",
]
