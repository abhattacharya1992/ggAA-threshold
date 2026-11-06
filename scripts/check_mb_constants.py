#!/usr/bin/env python3
"""Recompute the thirteen endpoint constants and their reductions.

Each MB[p,q,a,w] is evaluated by quadrature, its integer relation with the basis of
weight p + q is found with PSLQ, and the result is compared with the table in
ggaa/constants.py.  Writes data/constants/mb_constants.csv.
"""

import argparse
import json
import csv
import pathlib
import sys
from fractions import Fraction

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import mpmath as mp
import sympy as sp

from ggaa import constants
from ggaa.paths import CONSTANTS, RESULTS

SYMBOLS = {"zeta3": constants.ZETA3, "pi^2 ln2": constants.PI**2 * constants.LN2,
           "ln2^3": constants.LN2**3, "pi^4": constants.PI**4,
           "pi^2 ln2^2": constants.PI**2 * constants.LN2**2, "ln2^4": constants.LN2**4,
           "Li4(1/2)": constants.LI4_HALF, "zeta3 ln2": constants.ZETA3 * constants.LN2}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits", type=int, default=110)
    mp.mp.dps = parser.parse_args().digits

    CONSTANTS.mkdir(parents=True, exist_ok=True)
    agree = True
    residuals = []
    with open(CONSTANTS / "mb_constants.csv", "w", newline="") as handle:
        out = csv.writer(handle, lineterminator="\n")
        out.writerow(["p", "q", "a", "w", "value", "reduction"])
        for key, table_value in constants.REDUCTIONS.items():
            value, relation, residual = constants.identify(*key)
            found = sum(sp.Rational(c.numerator, c.denominator) * SYMBOLS[n] for n, c in relation.items())
            same = sp.simplify(found - table_value) == 0
            agree &= same
            label = "MB[{},{},{},{}]".format(*[str(Fraction(x)) for x in key])
            print(f"{label:16s} = {sp.sstr(table_value):70s} residual {mp.nstr(residual, 2):>8}  "
                  f"{'matches table' if same else 'DIFFERS FROM TABLE'}")
            out.writerow([str(Fraction(x)) for x in key] + [mp.nstr(value, 60), sp.sstr(table_value)])
            residuals.append(float(residual))
    print("\nAll reductions confirmed:", agree)
    json.dump({"digits": mp.mp.dps, "count": len(residuals), "max_residual": max(residuals), "confirmed": agree},
              open(RESULTS / "mb_constants.json", "w"), indent=1)
    if not agree:
        sys.exit(1)


if __name__ == "__main__":
    main()
