"""
Deterministic DNA sequence operations for BioSeqInsight.

These functions clean user input, then compute standard textbook statistics.
They do not discover motifs statistically (that is not MEME) and they do not
predict structure.
"""

from __future__ import annotations

import re
from typing import Iterable

from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt


DNA_BASES = set("ACGT")
DNA_AMBIGUOUS = set("ACGTN")
PROTEIN_LETTERS = set("ACDEFGHIKLMNPQRSTVWY")
STOP_SYMBOL = "*"

_STANDARD_TABLE = CodonTable.unambiguous_dna_by_id[1]


def strip_fasta_and_whitespace(text: str) -> str:
    """Remove FASTA headers, digits, and whitespace; keep letters only."""
    if not text:
        return ""
    lines = []
    for line in str(text).splitlines():
        if line.startswith(">"):
            continue
        lines.append(line)
    joined = "".join(lines)
    return "".join(ch for ch in joined if ch.isalpha())


def clean_dna(text: str) -> str:
    """Uppercase DNA, dropping characters outside A/C/G/T/N."""
    raw = strip_fasta_and_whitespace(text).upper().replace("U", "T")
    return "".join(ch for ch in raw if ch in DNA_AMBIGUOUS)


def clean_protein(text: str) -> str:
    """Uppercase protein, dropping non-standard residue letters (keeps X)."""
    raw = strip_fasta_and_whitespace(text).upper()
    return "".join(ch for ch in raw if ch in PROTEIN_LETTERS or ch == "X")


def looks_like_dna(text: str) -> bool:
    seq = strip_fasta_and_whitespace(text).upper()
    if len(seq) < 4:
        return False
    dna_frac = sum(ch in DNA_AMBIGUOUS for ch in seq) / len(seq)
    return dna_frac >= 0.85


def looks_like_protein(text: str) -> bool:
    seq = strip_fasta_and_whitespace(text).upper()
    if not seq:
        return False
    prot_frac = sum(ch in PROTEIN_LETTERS or ch == "X" for ch in seq) / len(seq)
    return prot_frac >= 0.85 and not looks_like_dna(text)


def gc_content(sequence: str) -> float:
    """GC% using A/C/G/T only (N excluded from numerator and denominator)."""
    seq = clean_dna(sequence)
    counted = "".join(ch for ch in seq if ch in DNA_BASES)
    if not counted:
        return 0.0
    return 100.0 * (counted.count("G") + counted.count("C")) / len(counted)


def nucleotide_frequency(sequence: str) -> dict[str, int]:
    seq = clean_dna(sequence)
    return {base: seq.count(base) for base in ("A", "T", "G", "C", "N")}


def reverse_sequence(sequence: str) -> str:
    return clean_dna(sequence)[::-1]


def reverse_complement(sequence: str) -> str:
    seq = clean_dna(sequence)
    if not seq:
        return ""
    return str(Seq(seq).reverse_complement())


def find_motif(sequence: str, motif: str) -> list[int]:
    """
    Exact, case-insensitive motif search.

    Returns 0-based start indices. This is string matching, not motif discovery.
    """
    seq = clean_dna(sequence)
    pat = clean_dna(motif)
    if not seq or not pat:
        return []
    positions = []
    start = 0
    while True:
        pos = seq.find(pat, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    return positions


def melting_temperature(sequence: str) -> dict:
    """
    Wallace rule for sequences <= 20 nt.
    Biopython nearest-neighbor (SantaLucia) for longer sequences.
    """
    seq = "".join(ch for ch in clean_dna(sequence) if ch in DNA_BASES)
    if not seq:
        return {"error": "No valid A/C/G/T bases found."}

    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")
    wallace = 2 * at + 4 * gc
    result = {
        "length": len(seq),
        "wallace_c": float(wallace),
        "method": "wallace",
        "tm_c": float(wallace),
    }
    if len(seq) <= 20:
        return result

    try:
        tm_nn = float(mt.Tm_NN(seq))
        result.update(
            {
                "method": "nearest_neighbor",
                "tm_c": tm_nn,
                "note": "Wallace is not meaningful for gene-length DNA.",
            }
        )
    except Exception as exc:
        result["error"] = f"Nearest-neighbor Tm failed: {exc}"
    return result


def dna_to_mrna(sequence: str) -> str:
    return clean_dna(sequence).replace("T", "U")


def translate_dna(sequence: str, table: int = 1) -> list[str]:
    """
    Translate frame 1 (nucleotides 0,1,2,...) with the NCBI standard table.

    Stop codons split the result into separate peptides. Trailing incomplete
    codons are ignored.
    """
    seq = "".join(ch for ch in clean_dna(sequence) if ch in DNA_BASES)
    if len(seq) < 3:
        return []
    usable = seq[: len(seq) - (len(seq) % 3)]
    protein = str(Seq(usable).translate(table=table, to_stop=False, cds=False))
    parts = [p for p in protein.split(STOP_SYMBOL) if p]
    return parts


def _codons(seq: str) -> Iterable[str]:
    for i in range(0, len(seq) - 2, 3):
        yield i, seq[i : i + 3]


def find_orfs(sequence: str, min_aa: int = 10, include_reverse: bool = True) -> list[dict]:
    """
    Find ATG-to-stop ORFs in the three forward frames, and optionally the
    three reverse frames. Coordinates are 0-based on the strand that was scanned.
    """
    fwd = "".join(ch for ch in clean_dna(sequence) if ch in DNA_BASES)
    strands = [("+", fwd)]
    if include_reverse:
        strands.append(("-", reverse_complement(fwd)))

    start_codons = set(_STANDARD_TABLE.start_codons) | {"ATG"}
    stop_codons = set(_STANDARD_TABLE.stop_codons)
    found: list[dict] = []

    for strand, seq in strands:
        for frame in range(3):
            in_orf = False
            aa: list[str] = []
            start_nt = None
            for nt_index, codon in _codons(seq[frame:]):
                abs_start = frame + nt_index
                if not in_orf:
                    if codon in start_codons:
                        in_orf = True
                        start_nt = abs_start
                        aa = ["M"]
                    continue
                if codon in stop_codons:
                    if len(aa) >= min_aa:
                        found.append(
                            {
                                "strand": strand,
                                "frame": frame + 1,
                                "start": start_nt,
                                "end": abs_start + 3,
                                "aa_length": len(aa),
                                "protein": "".join(aa),
                            }
                        )
                    in_orf = False
                    aa = []
                    start_nt = None
                else:
                    aa.append(_STANDARD_TABLE.forward_table.get(codon, "X"))
            # unterminated ORF at contig end is reported if long enough
            if in_orf and len(aa) >= min_aa:
                found.append(
                    {
                        "strand": strand,
                        "frame": frame + 1,
                        "start": start_nt,
                        "end": len(seq),
                        "aa_length": len(aa),
                        "protein": "".join(aa),
                        "unterminated": True,
                    }
                )
    return found


def format_gc_report(sequence: str) -> str:
    seq = clean_dna(sequence)
    counted = "".join(ch for ch in seq if ch in DNA_BASES)
    if not counted:
        return "No valid A/C/G/T bases found."
    freq = nucleotide_frequency(sequence)
    return (
        f"GC content: {gc_content(sequence):.2f}%\n"
        f"Length (A/C/G/T/N): {len(seq)} nt\n"
        f"Counts: A={freq['A']}, T={freq['T']}, G={freq['G']}, C={freq['C']}, N={freq['N']}"
    )


def format_frequency_report(sequence: str) -> str:
    seq = clean_dna(sequence)
    if not seq:
        return "No valid DNA bases found."
    freq = nucleotide_frequency(sequence)
    n = max(len(seq), 1)
    lines = [f"Length: {len(seq)} nt"]
    for base in ("A", "T", "G", "C", "N"):
        pct = 100.0 * freq[base] / n
        lines.append(f"{base}: {freq[base]} ({pct:.2f}%)")
    return "\n".join(lines)


def format_motif_report(sequence: str, motif: str) -> str:
    seq = clean_dna(sequence)
    pat = clean_dna(motif)
    if not pat:
        return "Enter a motif (A/C/G/T)."
    if not seq:
        return "No valid DNA sequence."
    hits = find_motif(seq, pat)
    if not hits:
        return f"Motif {pat} not found in {len(seq)} nt."
    one_based = [i + 1 for i in hits]
    return (
        f"Motif {pat}: {len(hits)} hit(s) in {len(seq)} nt\n"
        f"0-based starts: {hits}\n"
        f"1-based starts: {one_based}\n"
        "Note: exact string search, not MEME / PWM motif discovery."
    )


def format_tm_report(sequence: str) -> str:
    data = melting_temperature(sequence)
    if "error" in data and "tm_c" not in data:
        return data["error"]
    if data.get("method") == "wallace":
        return (
            f"Tm: {data['tm_c']:.1f} °C\n"
            f"Method: Wallace rule  Tm = 2(A+T)+4(G+C)\n"
            f"Length: {data['length']} nt (Wallace is intended for ~14–20 nt primers)"
        )
    lines = [
        f"Tm: {data['tm_c']:.2f} °C",
        "Method: nearest-neighbor (SantaLucia implementation in Biopython)",
        f"Length: {data['length']} nt",
        f"Wallace value on this sequence would be {data['wallace_c']:.0f} °C and is not usable.",
    ]
    if "error" in data:
        lines.append(data["error"])
    return "\n".join(lines)


def format_translate_report(sequence: str) -> str:
    seq = clean_dna(sequence)
    if looks_like_protein(sequence) and not looks_like_dna(sequence):
        return "This looks like a protein sequence. Use the 3D Structure tab."
    parts = translate_dna(seq)
    if not parts:
        return "Nothing to translate (need at least one complete codon of A/C/G/T)."
    blocks = [f"Frame +1 translation ({len(parts)} peptide(s); stops split peptides):\n"]
    for i, pep in enumerate(parts, 1):
        blocks.append(f">peptide_{i} len={len(pep)}\n{pep}")
    return "\n".join(blocks)


def format_orf_report(sequence: str, min_aa: int = 10) -> str:
    seq = clean_dna(sequence)
    if not seq:
        return "No valid DNA sequence."
    orfs = find_orfs(seq, min_aa=min_aa, include_reverse=True)
    if not orfs:
        return (
            f"No ORFs found with at least {min_aa} amino acids "
            "(ATG start, standard genetic code, 6 frames)."
        )
    orfs_sorted = sorted(orfs, key=lambda o: (-o["aa_length"], o["strand"], o["start"]))
    lines = [
        f"Found {len(orfs_sorted)} ORF(s) with min length {min_aa} aa.",
        "Coordinates are 0-based on the scanned strand.\n",
    ]
    for i, orf in enumerate(orfs_sorted, 1):
        flag = " (no stop codon)" if orf.get("unterminated") else ""
        lines.append(
            f"ORF {i}: strand {orf['strand']} frame {orf['frame']} "
            f"nt {orf['start']}-{orf['end']} ({orf['aa_length']} aa){flag}"
        )
        lines.append(orf["protein"])
        lines.append("")
    return "\n".join(lines).rstrip()
