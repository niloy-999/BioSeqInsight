"""Project (session) management.

A BioSeqInsight project is a self-describing directory that holds everything
needed to re-read or re-run an analysis:

::

    project/
      project.json      manifest: name, created/modified, software version
      metadata.json     the exact settings the analysis ran under
      sequences/        input FASTA files, copied in
      structures/       retrieved coordinate files
      results/          exported tables and JSON results
      logs/             the run log

The point is transferability. A project directory can be zipped and sent to a
collaborator, attached to a manuscript, or archived alongside a dataset, and
the recipient can see exactly which sequences were analysed, which structures
were retrieved from which resources, what identity each hit had and which
software version and settings produced the numbers.

Projects are plain files and plain JSON on purpose: no database, no pickled
Python objects, nothing that stops being readable when BioSeqInsight changes.
"""

from __future__ import annotations

import json
import shutil
import zipfile
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .. import __version__
from ..config.settings import Settings, get_settings
from ..core.alphabet import FastaRecord, parse_fasta, write_fasta
from ..io.export import to_json, write_csv, write_html_report, write_json
from ..models.results import BatchRow
from ..services.logging_setup import get_logger

logger = get_logger("project")

MANIFEST_NAME = "project.json"
METADATA_NAME = "metadata.json"
SUBDIRECTORIES = ("sequences", "structures", "results", "logs")
PROJECT_FORMAT_VERSION = 1


class ProjectError(RuntimeError):
    """Raised when a project cannot be created, opened or written."""


@dataclass
class ProjectManifest:
    """The contents of ``project.json``."""

    name: str
    created_at: str
    modified_at: str
    software_version: str
    format_version: int = PROJECT_FORMAT_VERSION
    description: str = ""
    sequence_files: list[str] = field(default_factory=list)
    structure_files: list[str] = field(default_factory=list)
    result_files: list[str] = field(default_factory=list)
    record_count: int = 0
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class Project:
    """An open project directory."""

    def __init__(self, path: str | Path, manifest: ProjectManifest, settings: Settings):
        self.path = Path(path)
        self.manifest = manifest
        self.settings = settings

    # -- lifecycle ----------------------------------------------------------

    @classmethod
    def create(
        cls,
        path: str | Path,
        name: str = "",
        description: str = "",
        settings: Settings | None = None,
        overwrite: bool = False,
    ) -> Project:
        """Create a new project directory."""
        target = Path(path)
        if target.exists():
            if not overwrite and any(target.iterdir()):
                raise ProjectError(
                    f"{target} already exists and is not empty. "
                    "Pass overwrite=True to reuse it."
                )
        target.mkdir(parents=True, exist_ok=True)
        for name_ in SUBDIRECTORIES:
            (target / name_).mkdir(exist_ok=True)

        now = _now()
        manifest = ProjectManifest(
            name=name or target.name,
            created_at=now,
            modified_at=now,
            software_version=__version__,
            description=description,
        )
        project = cls(target, manifest, settings or get_settings())
        project.save()
        logger.info("Created project %s", target)
        return project

    @classmethod
    def open(cls, path: str | Path, settings: Settings | None = None) -> Project:
        """Open an existing project directory."""
        target = Path(path)
        manifest_path = target / MANIFEST_NAME
        if not manifest_path.is_file():
            raise ProjectError(f"{target} is not a BioSeqInsight project ({MANIFEST_NAME} missing).")
        try:
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ProjectError(f"Cannot read {manifest_path}: {exc}") from exc

        format_version = int(data.get("format_version", 1))
        if format_version > PROJECT_FORMAT_VERSION:
            raise ProjectError(
                f"Project {target} uses format version {format_version}, but this "
                f"BioSeqInsight ({__version__}) understands up to {PROJECT_FORMAT_VERSION}. "
                "Upgrade BioSeqInsight to open it."
            )
        known = set(ProjectManifest.__dataclass_fields__)
        manifest = ProjectManifest(**{k: v for k, v in data.items() if k in known})
        loaded_settings = settings or _settings_from(target) or get_settings()
        return cls(target, manifest, loaded_settings)

    def save(self) -> None:
        """Write the manifest and the settings snapshot."""
        self.manifest.modified_at = _now()
        self.manifest.software_version = __version__
        try:
            (self.path / MANIFEST_NAME).write_text(
                json.dumps(self.manifest.to_dict(), indent=2), encoding="utf-8"
            )
            (self.path / METADATA_NAME).write_text(
                json.dumps(
                    {
                        "settings": self.settings.to_dict(),
                        "settings_source": self.settings.source,
                        "software_version": __version__,
                        "saved_at": self.manifest.modified_at,
                    },
                    indent=2,
                    sort_keys=True,
                ),
                encoding="utf-8",
            )
        except OSError as exc:
            raise ProjectError(f"Cannot write project files in {self.path}: {exc}") from exc

    # -- contents -----------------------------------------------------------

    def add_sequences(
        self, records: Sequence[FastaRecord], filename: str = "input.fasta"
    ) -> str:
        """Copy sequence records into the project."""
        target = self.path / "sequences" / filename
        target.write_text(write_fasta(records), encoding="utf-8")
        if filename not in self.manifest.sequence_files:
            self.manifest.sequence_files.append(filename)
        self.manifest.record_count = sum(
            len(parse_fasta((self.path / "sequences" / f).read_text(encoding="utf-8")))
            for f in self.manifest.sequence_files
            if (self.path / "sequences" / f).is_file()
        )
        self.save()
        return str(target)

    def import_fasta(self, source: str | Path) -> str:
        """Copy an existing FASTA file into the project unchanged."""
        source_path = Path(source)
        if not source_path.is_file():
            raise ProjectError(f"No such FASTA file: {source_path}")
        return self.add_sequences(
            parse_fasta(source_path.read_text(encoding="utf-8")), source_path.name
        )

    def read_sequences(self) -> list[FastaRecord]:
        """Every sequence record stored in the project."""
        records: list[FastaRecord] = []
        for filename in self.manifest.sequence_files:
            path = self.path / "sequences" / filename
            if path.is_file():
                records.extend(parse_fasta(path.read_text(encoding="utf-8")))
        return records

    def adopt_structure(self, source: str | Path) -> str | None:
        """Copy a retrieved coordinate file into the project."""
        source_path = Path(source)
        if not source_path.is_file():
            return None
        target = self.path / "structures" / source_path.name
        if source_path.resolve() != target.resolve():
            shutil.copy2(source_path, target)
        if source_path.name not in self.manifest.structure_files:
            self.manifest.structure_files.append(source_path.name)
        return str(target)

    def save_results(
        self,
        rows: Sequence[BatchRow],
        basename: str = "results",
        formats: Sequence[str] = ("csv", "json", "html"),
        adopt_structures: bool = True,
    ) -> dict[str, str]:
        """Write batch results into ``results/`` and adopt their structures."""
        written: dict[str, str] = {}
        results_dir = self.path / "results"
        if "csv" in formats:
            written["csv"] = write_csv(rows, results_dir / f"{basename}.csv")
        if "json" in formats:
            written["json"] = write_json(
                {
                    "software_version": __version__,
                    "generated_at": _now(),
                    "settings": self.settings.to_dict(),
                    "rows": [row.to_dict() for row in rows],
                },
                results_dir / f"{basename}.json",
            )
        if "html" in formats:
            written["html"] = write_html_report(
                rows,
                results_dir / f"{basename}.html",
                title=f"BioSeqInsight report - {self.manifest.name}",
            )

        if adopt_structures:
            for row in rows:
                for structure in row.structures:
                    if structure.download_path:
                        self.adopt_structure(structure.download_path)

        for path in written.values():
            filename = Path(path).name
            if filename not in self.manifest.result_files:
                self.manifest.result_files.append(filename)
        self.save()
        return written

    def write_note(self, text: str) -> None:
        self.manifest.notes.append(f"{_now()} {text}")
        self.save()

    # -- transfer -----------------------------------------------------------

    def export_zip(self, destination: str | Path | None = None) -> str:
        """Zip the whole project for transfer to another researcher."""
        target = Path(destination) if destination else self.path.with_suffix(".bsiproj.zip")
        target.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED) as archive:
            for item in sorted(self.path.rglob("*")):
                if item.is_file():
                    archive.write(item, item.relative_to(self.path.parent))
        logger.info("Exported project to %s", target)
        return str(target)

    @classmethod
    def import_zip(cls, archive_path: str | Path, destination: str | Path) -> Project:
        """Restore a project exported with :meth:`export_zip`."""
        destination = Path(destination)
        destination.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(archive_path) as archive:
            for member in archive.namelist():
                if member.startswith("/") or ".." in Path(member).parts:
                    raise ProjectError(f"Refusing to extract unsafe path {member!r}.")
            archive.extractall(destination)
        candidates = [destination] + [p.parent for p in destination.rglob(MANIFEST_NAME)]
        for candidate in candidates:
            if (candidate / MANIFEST_NAME).is_file():
                return cls.open(candidate)
        raise ProjectError(f"No {MANIFEST_NAME} found inside {archive_path}.")

    def summary(self) -> dict[str, Any]:
        return {
            "name": self.manifest.name,
            "path": str(self.path),
            "created_at": self.manifest.created_at,
            "modified_at": self.manifest.modified_at,
            "software_version": self.manifest.software_version,
            "records": self.manifest.record_count,
            "sequence_files": len(self.manifest.sequence_files),
            "structure_files": len(self.manifest.structure_files),
            "result_files": len(self.manifest.result_files),
        }

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"<Project {self.manifest.name!r} at {self.path}>"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _settings_from(path: Path) -> Settings | None:
    metadata = path / METADATA_NAME
    if not metadata.is_file():
        return None
    try:
        data = json.loads(metadata.read_text(encoding="utf-8"))
        return Settings.from_dict(data.get("settings", {}), source=str(metadata))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        logger.warning("Ignoring unreadable %s: %s", metadata, exc)
        return None


def list_projects(root: str | Path) -> list[dict[str, Any]]:
    """Summaries of every project directory under ``root``."""
    base = Path(root)
    if not base.is_dir():
        return []
    found = []
    for manifest in sorted(base.glob(f"*/{MANIFEST_NAME}")):
        try:
            found.append(Project.open(manifest.parent).summary())
        except ProjectError:  # pragma: no cover - skip corrupt directories
            continue
    return found


__all__ = [
    "MANIFEST_NAME",
    "PROJECT_FORMAT_VERSION",
    "Project",
    "ProjectError",
    "ProjectManifest",
    "list_projects",
    "to_json",
]
