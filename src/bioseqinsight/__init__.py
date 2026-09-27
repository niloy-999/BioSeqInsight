"""BioSeqInsight: an integrated, validated desktop workflow for DNA sequence
analysis and protein structure-resource retrieval.

The package is organised in layers:

``bioseqinsight.core``
    Deterministic, dependency-free sequence and protein computation.
``bioseqinsight.structures``
    Providers for external structural resources plus a unified manager,
    result model and identity-validation layer.
``bioseqinsight.services``
    Cross-cutting infrastructure: HTTP transport, retry policy, disk cache,
    logging.
``bioseqinsight.workflows``
    Batch analysis and project (session) management.
``bioseqinsight.io``
    Export to CSV/TSV/JSON/FASTA/HTML.
``bioseqinsight.gui``
    Tkinter front-end. The GUI only calls the public API above; it performs
    no computation of its own.

Nothing in ``core`` imports from ``gui``, ``structures`` or ``services``, so
every computational result in this package can be reproduced head-lessly and
tested without a network connection.
"""

from __future__ import annotations

__version__ = "2.0.2"
__author__ = "Tasnim Ul Islam, Md. Bayazid Hossen"
__license__ = "MIT"

__all__ = ["__author__", "__license__", "__version__"]
