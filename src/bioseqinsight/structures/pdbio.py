"""Reading PDB coordinate files.

Only what BioSeqInsight actually needs is implemented: chain sequences (so a
retrieved structure can be compared with the query), B-factor statistics (so
AlphaFold and ESMFold pLDDT can be reported), and the cross-references that
say which UniProt accession a structure claims to represent.

A full structural parser is out of scope. Users who need one already have
Biopython's ``Bio.PDB``, and BioSeqInsight does not attempt to replace it.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

THREE_TO_ONE: dict[str, str] = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    # Common modified residues mapped to their parent amino acid.
    "MSE": "M",
    "SEC": "U",
    "PYL": "O",
    "HYP": "P",
    "SEP": "S",
    "TPO": "T",
    "PTR": "Y",
    "CSO": "C",
    "CME": "C",
    "MLY": "K",
    "KCX": "K",
    "LLP": "K",
}

_DBREF = re.compile(r"^DBREF[12]?\s")
_UNIPROT_ACCESSION = re.compile(
    r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b"
)


def find_accessions(text: str) -> list[str]:
    """Every UniProt-shaped accession in a block of text, in order, deduplicated."""
    seen: set[str] = set()
    out: list[str] = []
    for match in _UNIPROT_ACCESSION.findall((text or "").upper()):
        if match not in seen:
            seen.add(match)
            out.append(match)
    return out


@dataclass
class Chain:
    """One polymer chain extracted from a coordinate file."""

    chain_id: str
    sequence: str
    residue_count: int
    mean_bfactor: float | None = None
    source: str = "atom"  # "atom" or "seqres"


@dataclass
class ParsedStructure:
    """The subset of a PDB file BioSeqInsight uses."""

    chains: list[Chain] = field(default_factory=list)
    mean_bfactor: float | None = None
    bfactor_range: str | None = None  # "0_1", "0_100" or "other"
    accessions: list[str] = field(default_factory=list)
    title: str = ""
    header_id: str = ""
    atom_count: int = 0
    has_coordinates: bool = False

    @property
    def longest_chain(self) -> Chain | None:
        return max(self.chains, key=lambda c: len(c.sequence), default=None)

    def best_matching_chain(self, query: str) -> Chain | None:
        """Chain whose length is closest to the query, preferring exact hits."""
        if not self.chains:
            return None
        target = len(query or "")
        exact = [c for c in self.chains if c.sequence == query]
        if exact:
            return exact[0]
        if not target:
            return self.longest_chain
        return min(self.chains, key=lambda c: abs(len(c.sequence) - target))


def parse_pdb(text: str) -> ParsedStructure:
    """Parse coordinates, chain sequences, B-factors and cross-references."""
    structure = ParsedStructure()
    if not text:
        return structure

    seqres: dict[str, list[str]] = {}
    atom_seq: dict[str, list[str]] = {}
    atom_bfactors: dict[str, list[float]] = {}
    seen_residues: set[tuple[str, str]] = set()
    all_bfactors: list[float] = []
    accessions: list[str] = []

    for line in text.splitlines():
        record = line[:6]
        if record == "HEADER":
            structure.header_id = line[62:66].strip()
        elif record == "TITLE ":
            structure.title = (structure.title + " " + line[10:].strip()).strip()
        elif _DBREF.match(line) and "UNP" in line:
            accessions.extend(find_accessions(line))
        elif record == "SEQRES":
            chain_id = line[11:12].strip() or "A"
            residues = line[19:].split()
            seqres.setdefault(chain_id, []).extend(residues)
        elif record in ("ATOM  ", "HETATM"):
            structure.atom_count += 1
            structure.has_coordinates = True
            altloc = line[16:17].strip()
            if altloc not in ("", "A"):
                continue
            resname = line[17:20].strip().upper()
            chain_id = line[21:22].strip() or "A"
            res_seq = line[22:27].strip()
            atom_name = line[12:16].strip()
            try:
                bfactor = float(line[60:66])
            except (ValueError, IndexError):
                bfactor = None
            if record == "ATOM  " and bfactor is not None:
                all_bfactors.append(bfactor)
            if atom_name == "CA" and resname in THREE_TO_ONE:
                key = (chain_id, res_seq)
                if key in seen_residues:
                    continue
                seen_residues.add(key)
                atom_seq.setdefault(chain_id, []).append(THREE_TO_ONE[resname])
                if bfactor is not None:
                    atom_bfactors.setdefault(chain_id, []).append(bfactor)

    chain_ids = sorted(set(atom_seq) | set(seqres))
    for chain_id in chain_ids:
        from_atoms = "".join(atom_seq.get(chain_id, []))
        from_seqres = "".join(THREE_TO_ONE.get(res, "X") for res in seqres.get(chain_id, []))
        # Prefer SEQRES when it is present and at least as complete: it
        # includes residues that are disordered and therefore absent from the
        # coordinates, which matters when comparing against a full-length
        # UniProt sequence.
        use_seqres = len(from_seqres) >= len(from_atoms) and bool(from_seqres)
        sequence = from_seqres if use_seqres else from_atoms
        if not sequence:
            continue
        bfs = atom_bfactors.get(chain_id, [])
        structure.chains.append(
            Chain(
                chain_id=chain_id,
                sequence=sequence,
                residue_count=len(sequence),
                mean_bfactor=(sum(bfs) / len(bfs)) if bfs else None,
                source="seqres" if use_seqres else "atom",
            )
        )

    if all_bfactors:
        mean = sum(all_bfactors) / len(all_bfactors)
        structure.mean_bfactor = mean
        structure.bfactor_range = _bfactor_range(all_bfactors)

    seen: set[str] = set()
    structure.accessions = [a for a in accessions if not (a in seen or seen.add(a))]
    return structure


def _bfactor_range(values: list[float]) -> str:
    """Report the numeric range of the B-factor column, nothing more.

    The column is overloaded: AlphaFold DB and ESM Atlas write per-atom pLDDT
    into it, while experimental entries carry real temperature factors, and
    both commonly fall in 0-100. The range alone therefore cannot tell them
    apart, so this function does not try. The *provider* knows what it asked
    for and is the one that labels the number as pLDDT or as a B-factor; see
    :func:`mean_confidence`.
    """
    maximum, minimum = max(values), min(values)
    if maximum <= 1.5 and minimum >= 0.0:
        return "0_1"
    if minimum >= 0.0 and maximum <= 100.0:
        return "0_100"
    return "other"


def mean_confidence(text: str, is_predicted: bool) -> tuple[float | None, str | None]:
    """Mean B-factor column value, labelled according to the source.

    ``is_predicted=True`` (AlphaFold DB, ESM Atlas) labels the value
    ``"plddt"`` and rescales a 0-1 column to 0-100. ``is_predicted=False``
    (experimental coordinates) labels it ``"bfactor"`` and never calls it a
    confidence score.
    """
    parsed = parse_pdb(text)
    if parsed.mean_bfactor is None:
        return None, None
    if not is_predicted:
        return parsed.mean_bfactor, "bfactor"
    if parsed.bfactor_range == "0_1":
        return parsed.mean_bfactor * 100.0, "plddt"
    return parsed.mean_bfactor, "plddt"


def looks_like_pdb(text: str) -> bool:
    """True when the payload contains coordinate records."""
    if not text:
        return False
    head = text[:200000]
    return "\nATOM  " in head or head.startswith("ATOM  ") or "\nHETATM" in head


__all__ = [
    "Chain",
    "ParsedStructure",
    "THREE_TO_ONE",
    "find_accessions",
    "looks_like_pdb",
    "mean_confidence",
    "parse_pdb",
]
