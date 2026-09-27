"""Structure visualisation helpers."""

from __future__ import annotations

import json
from pathlib import Path

from bioseqinsight.structures.visualization import (
    build_viewer_html,
    find_external_viewer,
    open_in_external_viewer,
    write_viewer_html,
)

MINIMAL = "ATOM      1  N   MET A   1      10.000  10.000  10.000  1.00 90.00           N\nEND\n"


class TestViewerPage:
    def test_page_is_well_formed(self):
        html = build_viewer_html(MINIMAL)
        assert html.startswith("<!DOCTYPE html>")
        assert "3Dmol" in html

    def test_coordinates_are_embedded_as_json(self):
        # v1.0 interpolated the PDB into a JavaScript template literal, so a
        # backtick or backslash in a REMARK produced a blank viewer. JSON
        # escaping cannot break that way.
        awkward = MINIMAL + "REMARK a backtick ` and a backslash \\ and ${injection}\n"
        html = build_viewer_html(awkward)
        assert json.dumps(awkward) in html
        assert "const pdbData = " in html

    def test_metadata_is_escaped(self):
        html = build_viewer_html(MINIMAL, title="<script>alert(1)</script>")
        assert "<script>alert(1)</script>" not in html
        assert "&lt;script&gt;" in html

    def test_mapping_level_is_displayed(self):
        html = build_viewer_html(MINIMAL, mapping="M4 Exact accession and sequence match")
        assert "M4 Exact" in html

    def test_written_to_a_given_path(self, tmp_path):
        target = tmp_path / "nested" / "view.html"
        path = write_viewer_html(MINIMAL, target)
        assert Path(path).is_file()

    def test_written_to_a_temporary_file_when_no_path_given(self):
        path = write_viewer_html(MINIMAL)
        assert Path(path).is_file() and path.endswith(".html")


class TestExternalViewer:
    def test_missing_file_is_reported(self, tmp_path):
        message = open_in_external_viewer(str(tmp_path / "absent.pdb"))
        assert "No coordinate file" in message

    def test_absent_preferred_viewer_falls_back_gracefully(self, tmp_path):
        pdb = tmp_path / "s.pdb"
        pdb.write_text(MINIMAL)
        message = open_in_external_viewer(str(pdb), preferred="definitely-not-installed")
        # One of three things happened, none of them an exception: a real
        # viewer was found; the OS's default-application handler took it
        # (macOS "open", Windows os.startfile -- neither names the viewer,
        # both name the file or say it was handed off); or the user was told
        # where the file is.
        assert any(marker in message for marker in (str(pdb), "Opened", "Handed", "default"))

    def test_find_external_viewer_returns_none_or_a_path(self):
        found = find_external_viewer()
        assert found is None or Path(found).exists()
