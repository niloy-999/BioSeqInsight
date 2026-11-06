"""Structure visualisation.

Visualisation is kept strictly separate from retrieval so that a coordinate
file is still useful when no viewer is installed, and so that no viewer is a
mandatory dependency of the sequence-analysis functionality.

Two routes are supported:

* a self-contained HTML page using 3Dmol.js, opened in the default browser;
* an external viewer (PyMOL, ChimeraX, or a user-configured command).

The HTML page embeds the coordinates as a JSON string literal rather than
interpolating them into a JavaScript template literal, which is what v1.0
did. Backticks or backslashes inside a PDB REMARK could break that template
and produce a blank viewer; JSON escaping cannot.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import webbrowser
from pathlib import Path

from ..services.logging_setup import get_logger

logger = get_logger("visualization")

VIEWER_CANDIDATES = ("pymol", "PyMOL", "chimerax", "ChimeraX", "vmd", "pymol.bat")

_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
<style>
  :root {{ color-scheme: dark; }}
  html, body {{ margin: 0; height: 100%; background: #0d1117; color: #e6edf3;
    font-family: ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif; }}
  #bar {{ padding: 8px 14px; font-size: 13px; border-bottom: 1px solid #30363d;
    display: flex; gap: 16px; flex-wrap: wrap; align-items: center; }}
  #bar b {{ font-weight: 600; }}
  #viewer {{ position: absolute; top: 37px; bottom: 0; left: 0; right: 0; }}
  .tag {{ background: #21262d; border-radius: 4px; padding: 2px 8px; }}
</style>
</head>
<body>
<div id="bar">
  <b>{title}</b>
  <span class="tag">{source}</span>
  <span class="tag">{confidence}</span>
  <span class="tag">{mapping}</span>
</div>
<div id="viewer"></div>
<script>
  const pdbData = {pdb_json};
  const viewer = $3Dmol.createViewer("viewer", {{ backgroundColor: "#0d1117" }});
  viewer.addModel(pdbData, "pdb");
  viewer.setStyle({{}}, {{ cartoon: {{ color: "spectrum" }} }});
  viewer.zoomTo();
  viewer.render();
</script>
</body>
</html>
"""


def build_viewer_html(
    pdb_text: str,
    title: str = "BioSeqInsight structure",
    source: str = "",
    confidence: str = "",
    mapping: str = "",
) -> str:
    """Return a self-contained HTML page that renders the coordinates."""
    return _TEMPLATE.format(
        title=_escape(title),
        source=_escape(source or "local file"),
        confidence=_escape(confidence or "confidence not reported"),
        mapping=_escape(mapping or "identity not checked"),
        pdb_json=json.dumps(pdb_text),
    )


def _escape(text: str) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write_viewer_html(
    pdb_text: str, destination: str | os.PathLike[str] | None = None, **kwargs
) -> str:
    """Write the viewer page to disk and return its path."""
    html = build_viewer_html(pdb_text, **kwargs)
    if destination is None:
        handle = tempfile.NamedTemporaryFile("w", suffix=".html", delete=False, encoding="utf-8")
        handle.write(html)
        handle.close()
        return handle.name
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(html, encoding="utf-8")
    return str(path)


def open_in_browser(html_path: str) -> bool:
    """Open a viewer page; returns False when no browser could be launched."""
    try:
        return webbrowser.open(Path(html_path).absolute().as_uri())
    except Exception as exc:  # pragma: no cover - platform dependent
        logger.warning("Could not open a browser: %s", exc)
        return False


def find_external_viewer(preferred: str = "") -> str | None:
    """Locate a molecular viewer on PATH, honouring a configured preference."""
    if preferred:
        resolved = shutil.which(preferred)
        if resolved:
            return resolved
        logger.warning("Configured viewer %r was not found on PATH.", preferred)
    for candidate in VIEWER_CANDIDATES:
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
    return None


def open_in_external_viewer(pdb_path: str, preferred: str = "") -> str:
    """Open a coordinate file in an external viewer; returns a status message."""
    path = Path(pdb_path)
    if not path.is_file():
        return f"No coordinate file at {path}."

    viewer = find_external_viewer(preferred)
    if viewer:
        try:
            subprocess.Popen([viewer, str(path)])  # noqa: S603 - user-selected viewer
            return f"Opened {path.name} in {Path(viewer).name}."
        except OSError as exc:
            logger.warning("Launching %s failed: %s", viewer, exc)

    if sys.platform == "darwin":
        try:
            subprocess.run(["open", str(path)], check=False)  # noqa: S603, S607
            return f"Handed {path.name} to the default macOS application."
        except OSError:  # pragma: no cover - platform specific
            pass

    return (
        f"No molecular viewer was found on PATH. The coordinate file is at {path}; "
        "open it with PyMOL, ChimeraX or any viewer, or use the built-in 3Dmol.js view."
    )


__all__ = [
    "build_viewer_html",
    "find_external_viewer",
    "open_in_browser",
    "open_in_external_viewer",
    "write_viewer_html",
]
