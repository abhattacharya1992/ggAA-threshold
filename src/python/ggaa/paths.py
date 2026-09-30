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
FIGURES = ROOT / "manuscript" / "figures"
TABLES = ROOT / "manuscript" / "tables"
