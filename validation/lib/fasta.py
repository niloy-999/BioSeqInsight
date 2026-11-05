"""Minimal FASTA I/O for the validation package (no dependency on bioseqinsight)."""

from __future__ import annotations

from pathlib import Path


def write_fasta(path: str | Path, records: list[tuple[str, str, str]], width: int = 60) -> Path:
    """Write ``[(id, description, sequence), ...]`` to a wrapped FASTA file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for identifier, description, sequence in records:
            header = f">{identifier} {description}".strip() if description else f">{identifier}"
            fh.write(header + "\n")
            for i in range(0, len(sequence), width):
                fh.write(sequence[i:i + width] + "\n")
    return path


def read_fasta(path: str | Path) -> list[tuple[str, str, str]]:
    """Parse a FASTA file into ``[(id, description, sequence), ...]``."""
    path = Path(path)
    records: list[tuple[str, str, str]] = []
    identifier = description = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if identifier is not None:
                records.append((identifier, description or "", "".join(chunks)))
            header = line[1:].strip()
            parts = header.split(None, 1)
            identifier = parts[0] if parts else ""
            description = parts[1] if len(parts) > 1 else ""
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if identifier is not None:
        records.append((identifier, description or "", "".join(chunks)))
    return records


__all__ = ["write_fasta", "read_fasta"]
