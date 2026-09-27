"""Tkinter front end.

Importing this package does not import Tkinter: on a headless machine the
rest of BioSeqInsight must still work. Use :func:`launch` (which imports the
window lazily) to start the GUI.
"""

from __future__ import annotations


def launch(settings=None) -> int:
    """Start the graphical interface, importing Tkinter only when asked to."""
    try:
        from .main_window import launch as _launch
    except ImportError as exc:  # pragma: no cover - depends on the Python build
        print(
            "The graphical interface needs Tkinter, which is not available in this "
            f"Python installation ({exc}).\n"
            "On Debian/Ubuntu: sudo apt install python3-tk\n"
            "Everything else is available from the command line: "
            "python -m bioseqinsight --help"
        )
        return 1
    return _launch(settings)


__all__ = ["launch"]
