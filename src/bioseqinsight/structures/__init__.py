"""Structural-resource subsystem: providers, validation and the manager."""

from .alphafold import AlphaFoldProvider
from .base import StructureProvider, classify_query, is_pdb_id, is_uniprot_accession
from .esmatlas import ESMAtlasProvider
from .manager import StructureManager, StructureQuery, StructureSearchOutcome
from .pdbio import parse_pdb
from .rcsb import RCSBProvider
from .uniprot import UniProtClient
from .validation import classify, validate_structure

__all__ = [
    "AlphaFoldProvider",
    "ESMAtlasProvider",
    "RCSBProvider",
    "StructureManager",
    "StructureProvider",
    "StructureQuery",
    "StructureSearchOutcome",
    "UniProtClient",
    "classify",
    "classify_query",
    "is_pdb_id",
    "is_uniprot_accession",
    "parse_pdb",
    "validate_structure",
]
