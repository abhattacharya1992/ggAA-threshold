"""Parse the extracted Nf^2 pieces into SymPy expressions."""

import pathlib
import pickle
import re

from sympy.parsing.mathematica import parse_mathematica


ROOT = pathlib.Path(__file__).parent

for label in ("P1", "P2"):
    text = (ROOT / f"{label}_nf2.wl").read_text()

    # The source uses N for the colour factor.  Rename it before passing the
    # expression to SymPy, where N otherwise has an inconvenient meaning.
    expression = parse_mathematica(re.sub(r"\bN\b", "nc", text))

    output = ROOT / f"{label}_nf2_expr.pkl"
    output.write_bytes(pickle.dumps(expression))
    print("Parsed", label)
