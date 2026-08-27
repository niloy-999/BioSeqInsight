"""
Protein-level calculations and ESM Atlas API structure requests.

Important: predict_structure() calls a remote service. It does not load
ESMFold weights on this computer.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import webbrowser

import requests
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from sequence_operations import clean_protein, looks_like_dna


ESMATLAS_FOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
DEFAULT_PDB_NAME = "predicted_structure.pdb"
STRUCTURE_DIR = "structures"
LAST_OPENED_PDB = None
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{uid}"
ALPHAFOLD_FILE = "https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_{ver}.pdb"
RCSB_PDB = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
PDBE_UNIPROT_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
UNIPROT_JSON_URL = "https://rest.uniprot.org/uniprotkb/{uid}.json"

UNIPROT_RE = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$",
    re.I,
)
PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")

# Kyte–Doolittle hydropathy (Kyte & Doolittle, 1982)
KYTE_DOOLITTLE = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

# Simplified single-residue helix/sheet/coil propensities.
# These are teaching approximations, not full Chou–Fasman or GOR.
HELIX_P = {
    "A": 1.42,
    "C": 0.70,
    "D": 1.01,
    "E": 1.51,
    "F": 1.13,
    "G": 0.57,
    "H": 1.00,
    "I": 1.08,
    "K": 1.16,
    "L": 1.21,
    "M": 1.45,
    "N": 0.67,
    "P": 0.57,
    "Q": 1.11,
    "R": 0.98,
    "S": 0.77,
    "T": 0.83,
    "V": 1.06,
    "W": 1.08,
    "Y": 0.69,
}
SHEET_P = {
    "A": 0.83,
    "C": 1.19,
    "D": 0.54,
    "E": 0.37,
    "F": 1.38,
    "G": 0.75,
    "H": 0.87,
    "I": 1.60,
    "K": 0.74,
    "L": 1.30,
    "M": 1.05,
    "N": 0.89,
    "P": 0.55,
    "Q": 1.10,
    "R": 0.93,
    "S": 0.75,
    "T": 1.19,
    "V": 1.70,
    "W": 1.37,
    "Y": 1.47,
}
COIL_P = {
    "A": 0.66,
    "C": 1.19,
    "D": 1.46,
    "E": 0.74,
    "F": 0.60,
    "G": 1.52,
    "H": 0.87,
    "I": 0.47,
    "K": 1.20,
    "L": 0.59,
    "M": 0.60,
    "N": 1.56,
    "P": 1.52,
    "Q": 0.92,
    "R": 0.93,
    "S": 1.43,
    "T": 1.20,
    "V": 0.61,
    "W": 0.60,
    "Y": 1.14,
}


def _require_protein(sequence: str) -> str:
    if looks_like_dna(sequence):
        raise ValueError(
            "This input looks like DNA. Translate it on the Sequence Analysis "
            "tab, then paste the protein here."
        )
    seq = clean_protein(sequence)
    if not seq:
        raise ValueError("No valid amino-acid letters found.")
    return seq


def calculate_hydrophobicity(sequence: str) -> float:
    """Mean Kyte–Doolittle hydropathy (GRAVY-like)."""
    seq = _require_protein(sequence)
    values = [KYTE_DOOLITTLE[aa] for aa in seq if aa in KYTE_DOOLITTLE]
    if not values:
        return 0.0
    return sum(values) / len(values)


def calculate_molecular_weight(sequence: str) -> float:
    """Average isotopic molecular weight in Daltons (peptide, water loss included)."""
    seq = _require_protein(sequence)
    # ProteinAnalysis ignores X; drop unknowns for a defined mass
    usable = "".join(ch for ch in seq if ch in KYTE_DOOLITTLE)
    if not usable:
        raise ValueError("Sequence has no standard amino acids.")
    return float(ProteinAnalysis(usable).molecular_weight())


def predict_secondary_structure(sequence: str, window: int = 5) -> dict:
    """
    Smoothed residue-wise propensity sketch.

    Not a full Chou–Fasman or GOR implementation. Use it only as a rough
    visual, not as a structural assignment for publication.
    """
    seq = _require_protein(sequence)
    raw = []
    for aa in seq:
        scores = {
            "H": HELIX_P.get(aa, 0.0),
            "E": SHEET_P.get(aa, 0.0),
            "C": COIL_P.get(aa, 0.0),
        }
        raw.append(max(scores, key=scores.get))

    half = max(window // 2, 0)
    smoothed = []
    for i in range(len(raw)):
        chunk = raw[max(0, i - half) : i + half + 1]
        smoothed.append(max(set(chunk), key=chunk.count))

    visual = "".join({"H": "█", "E": "≈", "C": "."}[s] for s in smoothed)
    n = len(smoothed) or 1
    helix = 100.0 * smoothed.count("H") / n
    sheet = 100.0 * smoothed.count("E") / n
    coil = 100.0 * smoothed.count("C") / n
    return {
        "length": len(seq),
        "visual": visual,
        "helix_pct": helix,
        "sheet_pct": sheet,
        "coil_pct": coil,
        "method": f"smoothed single-residue propensities (window={window})",
    }


def format_hydrophobicity_report(sequence: str) -> str:
    try:
        score = calculate_hydrophobicity(sequence)
    except ValueError as exc:
        return str(exc)
    seq = clean_protein(sequence)
    label = "hydrophobic" if score > 0 else "hydrophilic"
    return (
        f"Mean Kyte–Doolittle hydropathy: {score:.3f}\n"
        f"Interpretation: {label} overall (0 is a common midpoint)\n"
        f"Length used: {len(seq)} aa"
    )


def format_mw_report(sequence: str) -> str:
    try:
        mass = calculate_molecular_weight(sequence)
    except ValueError as exc:
        return str(exc)
    seq = clean_protein(sequence)
    return (
        f"Molecular weight: {mass:.2f} Da ({mass / 1000:.3f} kDa)\n"
        f"Length: {len(seq)} aa\n"
        "Computed with Biopython ProteinAnalysis (average isotopic mass)."
    )


def format_ss_report(sequence: str) -> str:
    try:
        data = predict_secondary_structure(sequence)
    except ValueError as exc:
        return str(exc)
    return (
        "Simplified secondary-structure sketch "
        "(not full Chou–Fasman / GOR; not for publication claims)\n"
        f"Method: {data['method']}\n"
        f"Length: {data['length']} aa\n"
        f"Helix █ {data['helix_pct']:.1f}%   "
        f"Sheet ≈ {data['sheet_pct']:.1f}%   "
        f"Coil . {data['coil_pct']:.1f}%\n\n"
        f"{data['visual']}"
    )


def _looks_like_predicted_model(pdb_text: str, hint: str = "") -> bool:
    blob = f"{hint}\n{pdb_text[:2500]}".upper()
    markers = (
        "ALPHAFOLD",
        "ESMFOLD",
        "ESM FOLD",
        "P-LDDT",
        "PLDDT",
        "AF-MODEL",
        "AF-P",
        "AF-Q",
        "AF-O",
        "MODEL CONFIDENCE",
    )
    return any(marker in blob for marker in markers)


def _format_ca_score(pdb_text: str, hint: str = "") -> str:
    """Label pLDDT for predicted models and B-factors for crystal structures."""
    mean = _mean_plddt_from_pdb(pdb_text)
    if mean is None:
        return "CA B-column: not found"
    if _looks_like_predicted_model(pdb_text, hint):
        return f"Mean CA pLDDT: {mean:.2f} / 100 (predicted-model confidence)"
    return (
        f"Mean CA B-factor: {mean:.2f} Å² "
        "(crystallographic mobility; this is not pLDDT)"
    )


def _structure_path(filename: str) -> str:
    os.makedirs(STRUCTURE_DIR, exist_ok=True)
    return os.path.abspath(os.path.join(STRUCTURE_DIR, filename))


def _resolved_pdb_path(pdb_path: str) -> str:
    if os.path.dirname(pdb_path):
        folder = os.path.dirname(os.path.abspath(pdb_path))
        os.makedirs(folder, exist_ok=True)
        return os.path.abspath(pdb_path)
    return _structure_path(os.path.basename(pdb_path))


def _mean_plddt_from_pdb(pdb_text: str) -> float | None:
    """Read CA B-column. Predicted models use this for pLDDT; crystals use B-factors."""
    values = []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom = line[12:16].strip()
        if atom not in {"CA", "C"}:
            # ESMFold puts pLDDT on every atom; prefer CA when present
            continue
        try:
            values.append(float(line[60:66]))
        except ValueError:
            continue
    if not values:
        # fall back to all ATOM B-factors
        for line in pdb_text.splitlines():
            if line.startswith("ATOM"):
                try:
                    values.append(float(line[60:66]))
                except ValueError:
                    continue
    if not values:
        return None
    mean = sum(values) / len(values)
    # Some pipelines store 0–1, others 0–100
    if mean <= 1.5:
        mean *= 100.0
    return mean


def create_3d_structure_html(pdb_data: str) -> str:
    safe_pdb = pdb_data.replace("\\", "\\\\").replace("`", "\\`")
    html_content = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8"/>
  <title>BioSeqInsight structure view</title>
  <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
  <style>
    html, body, #container {{ margin: 0; height: 100%; background: #111; }}
  </style>
</head>
<body>
  <div id="container"></div>
  <script>
    const viewer = $3Dmol.createViewer("container", {{ backgroundColor: "black" }});
    viewer.addModel(`{safe_pdb}`, "pdb");
    viewer.setStyle({{}}, {{ cartoon: {{ color: "spectrum" }} }});
    viewer.zoomTo();
    viewer.render();
  </script>
</body>
</html>
"""
    handle = tempfile.NamedTemporaryFile(delete=False, suffix=".html", mode="w", encoding="utf-8")
    handle.write(html_content)
    handle.close()
    return handle.name


def _save_and_view(pdb_string: str, pdb_path: str, open_browser: bool, header: str) -> str:
    global LAST_OPENED_PDB
    abs_pdb = _resolved_pdb_path(pdb_path)
    with open(abs_pdb, "w", encoding="utf-8") as handle:
        handle.write(pdb_string)
    LAST_OPENED_PDB = abs_pdb
    html_path = create_3d_structure_html(pdb_string)
    if open_browser:
        webbrowser.open("file://" + html_path)
    return (
        f"{header}\n"
        f"{_format_ca_score(pdb_string, header + ' ' + abs_pdb)}\n"
        f"PDB saved to: {abs_pdb}\n"
        f"Viewer HTML: {html_path}"
    )


def lookup_sequence_matches(sequence: str, identity_cutoff: float = 0.9, max_pdb: int = 8) -> dict:
    """
    Find PDB and UniProt IDs for a protein sequence via RCSB sequence search
    plus PDBe/UniProt metadata.
    """
    seq = _require_protein(sequence)
    payload = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": 1e-6,
                "identity_cutoff": identity_cutoff,
                "sequence_type": "protein",
                "value": seq,
            },
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": max_pdb},
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }
    try:
        response = requests.post(RCSB_SEARCH_URL, json=payload, timeout=45)
        response.raise_for_status()
        hits = response.json().get("result_set") or []
    except requests.exceptions.RequestException as exc:
        raise RuntimeError(f"RCSB sequence search failed: {exc}") from exc

    pdb_ids = [str(hit.get("identifier", "")).upper() for hit in hits if hit.get("identifier")]
    info = {
        "query_length": len(seq),
        "pdb_ids": pdb_ids,
        "uniprot": None,
        "uniprot_entry": None,
        "protein_name": None,
        "organism": None,
        "uniprot_length": None,
    }
    for pdb_id in pdb_ids[:4]:
        try:
            mapped = requests.get(PDBE_UNIPROT_URL.format(pdb_id=pdb_id.lower()), timeout=20)
            mapped.raise_for_status()
            block = mapped.json().get(pdb_id.lower(), {}).get("UniProt") or {}
        except requests.exceptions.RequestException:
            continue
        if not block:
            continue
        uid = sorted(block.keys(), key=lambda k: -len(block[k].get("mappings") or []))[0]
        rec = block[uid]
        info["uniprot"] = uid
        info["uniprot_entry"] = rec.get("name") or rec.get("identifier")
        break

    if info["uniprot"]:
        try:
            meta = requests.get(UNIPROT_JSON_URL.format(uid=info["uniprot"]), timeout=20)
            if meta.status_code == 200:
                js = meta.json()
                info["uniprot_entry"] = js.get("uniProtkbId") or info["uniprot_entry"]
                info["organism"] = (js.get("organism") or {}).get("scientificName")
                info["uniprot_length"] = (js.get("sequence") or {}).get("length")
                info["protein_name"] = (
                    ((js.get("proteinDescription") or {}).get("recommendedName") or {})
                    .get("fullName") or {}
                ).get("value")
        except requests.exceptions.RequestException:
            pass
    return info


def format_match_report(info: dict) -> str:
    if not info.get("pdb_ids") and not info.get("uniprot"):
        return "No close PDB / UniProt match was found for this sequence."
    lines = [f"Query length: {info.get('query_length')} aa"]
    if info.get("uniprot"):
        lines.append(f"UniProt: {info['uniprot']}")
        if info.get("uniprot_entry"):
            lines.append(f"Entry: {info['uniprot_entry']}")
        if info.get("protein_name"):
            lines.append(f"Protein: {info['protein_name']}")
        if info.get("organism"):
            lines.append(f"Organism: {info['organism']}")
        if info.get("uniprot_length"):
            lines.append(f"UniProt length: {info['uniprot_length']} aa")
        lines.append(f"AlphaFold DB: https://alphafold.ebi.ac.uk/entry/{info['uniprot']}")
    if info.get("pdb_ids"):
        lines.append("PDB IDs: " + ", ".join(info["pdb_ids"]))
        lines.append(f"RCSB top hit: https://www.rcsb.org/structure/{info['pdb_ids'][0]}")
    return "\n".join(lines)


def find_and_download_matches(
    sequence: str,
    pdb_path: str = DEFAULT_PDB_NAME,
    open_browser: bool = True,
) -> str:
    """Look up IDs from a sequence and download AlphaFold DB + top experimental PDB."""
    try:
        info = lookup_sequence_matches(sequence)
    except (ValueError, RuntimeError) as exc:
        return str(exc)

    report = ["Sequence ID lookup", format_match_report(info), ""]
    downloaded = []
    view_path = None

    if info.get("uniprot"):
        af_path = _structure_path(f"AF-{info['uniprot']}.pdb")
        result = fetch_alphafold_db(info["uniprot"], pdb_path=af_path, open_browser=False)
        report.append(result)
        if os.path.isfile(af_path) and "Downloaded AlphaFold" in result:
            downloaded.append(af_path)
            view_path = af_path
    else:
        report.append("No UniProt accession found, so AlphaFold DB was skipped.")

    if info.get("pdb_ids"):
        pdb_id = info["pdb_ids"][0]
        exp_path = _structure_path(f"{pdb_id}.pdb")
        result = fetch_rcsb_pdb(pdb_id, pdb_path=exp_path, open_browser=False)
        report.append(result)
        if os.path.isfile(exp_path) and "Downloaded experimental" in result:
            downloaded.append(exp_path)
            if view_path is None:
                view_path = exp_path
    else:
        report.append("No PDB hit found, so RCSB download was skipped.")

    if view_path and os.path.isfile(view_path):
        with open(view_path, encoding="utf-8", errors="replace") as handle:
            view = _save_and_view(
                handle.read(),
                view_path,
                open_browser,
                f"Opened {os.path.basename(view_path)} in the viewer. "
                "Other downloads were kept as separate files.",
            )
        report.append(view)

    if downloaded:
        report.append("Saved files:\n" + "\n".join(downloaded))
        report.append(
            "Files are not copied onto predicted_structure.pdb, so earlier "
            "ESM Atlas predictions are left unchanged."
        )
    return "\n\n".join(part for part in report if part)


def fetch_alphafold_db(uniprot_id: str, pdb_path: str = DEFAULT_PDB_NAME, open_browser: bool = True) -> str:
    """Download a precomputed AlphaFold DB model by UniProt accession."""
    raw = uniprot_id.strip()
    uid = raw.split()[0].split("|")[-1].upper()
    if uid.startswith("AF-"):
        uid = uid.split("-")[1]
    if not UNIPROT_RE.match(uid):
        try:
            info = lookup_sequence_matches(raw)
        except (ValueError, RuntimeError) as exc:
            return str(exc)
        if not info.get("uniprot"):
            return (
                "No UniProt accession found for this sequence.\n"
                + format_match_report(info)
            )
        uid = info["uniprot"]
        prefix = format_match_report(info) + "\n\n"
    else:
        prefix = ""
    if pdb_path == DEFAULT_PDB_NAME:
        pdb_path = _structure_path(f"AF-{uid}.pdb")
    errors = []
    for ver in ("v6", "v4"):
        url = ALPHAFOLD_FILE.format(uid=uid, ver=ver)
        try:
            response = requests.get(url, timeout=30)
        except requests.exceptions.RequestException as exc:
            errors.append(f"{url}: {exc}")
            continue
        if response.status_code == 200 and "ATOM" in response.text:
            return _save_and_view(
                response.text,
                pdb_path,
                open_browser,
                f"{prefix}Downloaded AlphaFold DB model {ver} for {uid}\nSource: {url}\n"
                "This is a stored model, not a live ESMFold run.",
            )
        errors.append(f"{url}: HTTP {response.status_code}")
    try:
        meta = requests.get(ALPHAFOLD_API.format(uid=uid), timeout=30)
        errors.append(f"API {meta.status_code}: {meta.text[:180]}")
    except requests.exceptions.RequestException as exc:
        errors.append(str(exc))
    return "AlphaFold DB lookup failed:\n" + "\n".join(errors)


def fetch_rcsb_pdb(pdb_id: str, pdb_path: str = DEFAULT_PDB_NAME, open_browser: bool = True) -> str:
    """Download an experimental PDB from RCSB."""
    raw = pdb_id.strip()
    pid = raw.split()[0].upper()
    prefix = ""
    if not PDB_ID_RE.match(pid):
        try:
            info = lookup_sequence_matches(raw)
        except (ValueError, RuntimeError) as exc:
            return str(exc)
        if not info.get("pdb_ids"):
            return "No PDB entry found for this sequence.\n" + format_match_report(info)
        pid = info["pdb_ids"][0]
        prefix = format_match_report(info) + "\n\n"
    if pdb_path == DEFAULT_PDB_NAME:
        pdb_path = _structure_path(f"{pid}.pdb")
    url = RCSB_PDB.format(pdb_id=pid)
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
    except requests.exceptions.RequestException as exc:
        return f"RCSB download failed: {exc}"
    if "ATOM" not in response.text and "HETATM" not in response.text:
        return f"RCSB returned no coordinates for {pid}."
    return _save_and_view(
        response.text,
        pdb_path,
        open_browser,
        f"{prefix}Downloaded experimental structure {pid} from RCSB\nSource: {url}",
    )


def load_pdb_file(file_path: str, open_browser: bool = True) -> str:
    if not file_path or not os.path.isfile(file_path):
        return f"PDB file not found: {file_path}"
    with open(file_path, encoding="utf-8", errors="replace") as handle:
        pdb_string = handle.read()
    if "ATOM" not in pdb_string and "HETATM" not in pdb_string:
        return "That file has no ATOM records."
    source = os.path.abspath(file_path)
    return _save_and_view(
        pdb_string,
        source,
        open_browser,
        f"Opened local PDB without copying it: {source}",
    )


def _post_esmatlas(sequence: str, timeout: int) -> tuple[str | None, str]:
    headers_list = (
        {"Content-Type": "application/x-www-form-urlencoded"},
        {"Content-Type": "text/plain"},
    )
    last = ""
    for headers in headers_list:
        try:
            response = requests.post(
                ESMATLAS_FOLD_URL,
                headers=headers,
                data=sequence.encode("utf-8"),
                timeout=timeout,
            )
        except requests.exceptions.RequestException as exc:
            last = str(exc)
            continue
        if response.status_code == 200 and "ATOM" in response.text:
            return response.text, f"HTTP 200 via {headers['Content-Type']}"
        last = f"HTTP {response.status_code} {response.reason}"
        snippet = response.text[:180].replace("\n", " ")
        if snippet:
            last += f" — {snippet}"
        if response.status_code not in (408, 429, 500, 502, 503, 504):
            break
    return None, last


def predict_structure(
    sequence: str,
    pdb_path: str = DEFAULT_PDB_NAME,
    timeout: int = 90,
    open_browser: bool = True,
) -> str:
    """
    Try the ESM Atlas fold API with retries.

    504/timeout means the public server is overloaded. The app still cannot
    run ESMFold locally. Use fetch_alphafold_db() or a local PDB instead.
    """
    raw = sequence.strip()
    token = raw.split()[0] if raw else ""
    if UNIPROT_RE.match(token):
        return fetch_alphafold_db(token, pdb_path=pdb_path, open_browser=open_browser)
    if PDB_ID_RE.match(token) and len(raw) <= 6:
        return fetch_rcsb_pdb(token, pdb_path=pdb_path, open_browser=open_browser)

    try:
        seq = _require_protein(sequence)
    except ValueError as exc:
        return str(exc)

    if len(seq) > 400:
        return (
            f"Sequence is {len(seq)} aa. The public ESM Atlas API is unreliable "
            "above ~400 aa and often returns 504.\n"
            "Paste a UniProt accession (e.g. P69905) to fetch AlphaFold DB, "
            "or a PDB ID (e.g. 1UBQ), or load examples/1UBQ.pdb."
        )

    attempts = []
    timeouts = (timeout, timeout + 60, timeout + 120)
    for i, seconds in enumerate(timeouts, 1):
        pdb_text, detail = _post_esmatlas(seq, seconds)
        attempts.append(f"attempt {i} ({seconds}s): {detail}")
        if pdb_text:
            if pdb_path == DEFAULT_PDB_NAME:
                stamp = time.strftime("%Y%m%d-%H%M%S")
                pdb_path = _structure_path(f"esmfold_{len(seq)}aa_{stamp}.pdb")
            return _save_and_view(
                pdb_text,
                pdb_path,
                open_browser,
                "ESMFold prediction from the ESM Atlas fold API "
                f"({ESMATLAS_FOLD_URL}).\n"
                f"Length: {len(seq)} aa. This is a remote call, not local ESMFold.",
            )
        if i < len(timeouts):
            time.sleep(2)

    return (
        "ESM Atlas public fold API failed.\n"
        + "\n".join(attempts)
        + "\n\n504 Gateway Timeout means their server gave up, not that your "
        "sequence is invalid. That endpoint is a shared demo service and is "
        "often down.\n\n"
        "What works right now:\n"
        "1. Paste a UniProt ID (example P69905) and click Fetch AlphaFold DB.\n"
        "2. Paste a PDB ID (example 1UBQ) and click Fetch RCSB structure.\n"
        "3. Load examples/1UBQ.pdb to demo the viewer offline.\n"
        "4. Retry Predict later with a protein under ~150 aa."
    )


def open_in_bioviewer(pdb_file: str | None = None) -> str:
    """Open a PDB in BioViewer on macOS, else fall back to any available viewer."""
    pdb_file_path = pdb_file or LAST_OPENED_PDB or DEFAULT_PDB_NAME
    pdb_file_path = os.path.abspath(pdb_file_path)
    if not os.path.isfile(pdb_file_path):
        return (
            f"No PDB file at {pdb_file_path}. "
            "Run a prediction or download first. Files are saved under structures/."
        )

    if sys.platform == "darwin":
        try:
            subprocess.run(["open", "-a", "BioViewer", pdb_file_path], check=True)
            return f"Opened in BioViewer: {pdb_file_path}"
        except (subprocess.CalledProcessError, FileNotFoundError):
            subprocess.run(["open", pdb_file_path], check=False)
            return (
                "BioViewer app not available. Opened the PDB with the "
                f"default macOS handler: {pdb_file_path}"
            )

    for candidate in ("pymol", "pymol3", "chimerax", "pymol.bat"):
        if shutil.which(candidate):
            subprocess.Popen([candidate, pdb_file_path])
            return f"Opened with {candidate}: {pdb_file_path}"

    webbrowser.open("file://" + pdb_file_path)
    return (
        f"No BioViewer / PyMOL / ChimeraX found. PDB path:\n{pdb_file_path}\n"
        "Open that file in any molecular viewer."
    )
