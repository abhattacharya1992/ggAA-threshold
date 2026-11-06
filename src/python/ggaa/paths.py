"""Locations inside the repository."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data"
CONSTRUCTION = DATA / "construction"
SERIES = DATA / "series"
VALIDATION = DATA / "validation"
CONSTANTS = DATA / "constants"
INPUTS = ROOT / "inputs" / "arXiv-1909.08993" / "anc"
RESULTS = ROOT / "results"


def _manuscript():
    """manuscript/ in the repository; in the arXiv tarball the package is anc/ and the
    manuscript (main.tex) is the directory above it."""
    local = ROOT / "manuscript"
    if local.is_dir():
        return local
    if (ROOT.parent / "main.tex").exists():
        return ROOT.parent
    return local


MANUSCRIPT = _manuscript()
FIGURES = MANUSCRIPT / "figures"
TABLES = MANUSCRIPT / "tables"
