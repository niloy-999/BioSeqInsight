import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sequence_operations import (
    clean_dna,
    find_motif,
    find_orfs,
    gc_content,
    looks_like_dna,
    looks_like_protein,
    reverse_complement,
    translate_dna,
)


def test_gc_content_simple():
    assert gc_content("ATGC") == 50.0
    assert gc_content(">hdr\natgc\n") == 50.0
    assert gc_content("gggg") == 100.0


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"


def test_motif_positions():
    assert find_motif("ATGATG", "ATG") == [0, 3]


def test_translate_does_not_replace_g_with_c():
    # ATG GGG TAA  -> M G
    peptides = translate_dna("ATGGGG TAA")
    assert peptides == ["MG"]


def test_translate_known_peptide():
    # ATG TCT TCT AAA GTT AAA TAA -> MSSKVK
    peptides = translate_dna("ATGTCTTCTAAAGTTAAATAA")
    assert peptides == ["MSSKVK"]


def test_orf_finds_mini_gene():
    orfs = find_orfs("CCCATGTCTTCTAAAGTTAAATAACCC", min_aa=4, include_reverse=False)
    proteins = [o["protein"] for o in orfs]
    assert "MSSKVK" in proteins


def test_protein_vs_dna_guess():
    assert looks_like_dna("ATGCATGCATGC")
    assert looks_like_protein("MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG")
    assert not looks_like_protein("ATGCATGCATGCATGC")


def test_clean_ignores_headers():
    assert clean_dna(">gene\nAtgC\n") == "ATGC"
