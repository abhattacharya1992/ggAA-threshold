#!/usr/bin/env python3
"""Compare the final series with the 1080 validation records.

Two things are checked:
  * the Python evaluation of data/series reproduces the series values stored in the
    validation table (which were computed in Mathematica from the unreduced series);
  * the truncation errors with respect to the unexpanded remainders, which give
    Table II and the error slopes of Fig. 1.
Writes results/accuracy.csv and manuscript/tables/accuracy.tex (Table II).
"""

import csv
import math
import json
import pathlib
import sys
from collections import defaultdict
from fractions import Fraction

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import mpmath as mp
import numpy as np

from ggaa.paths import RESULTS, TABLES
from ggaa.series import ThresholdSeries
from ggaa.validation import records, relative_error

mp.mp.dps = 50


def evaluate(series, rec, order):
    name, beta, c = rec["quantity"], rec["beta"], rec["costheta"]
    i, nf = int(name[1]), name.split("_nf")[1]
    kwargs = dict(beta=mp.mpf(beta.numerator) / beta.denominator,
                  costh=mp.mpf(c.numerator) / c.denominator,
                  mA=rec["mA"], muR=rec["muR"], Nc=rec["Nc"], order=order)
    if nf == "5":
        return series(i, nf=5, **kwargs)
    return series.component(i, int(nf), **kwargs)


REPRODUCTION_TOLERANCE = 1e-30   # the stored series values keep 40 digits


def main():
    series = ThresholdSeries()
    worst_reproduction = 0
    errors = defaultdict(lambda: defaultdict(float))       # (quantity, order) -> beta -> max error
    for count, rec in enumerate(records(), 1):
        for order in (4, 6):
            value = evaluate(series, rec, order)
            saved = rec[f"series{order}"]
            worst_reproduction = max(worst_reproduction, relative_error(value, saved))
            key = (rec["quantity"], order)
            errors[key][rec["beta"]] = max(errors[key][rec["beta"]], float(relative_error(value, rec["exact"])))
    print(f"{count} records.  Largest relative difference between this evaluation and the "
          f"stored series values: {mp.nstr(worst_reproduction, 3)}")

    betas = sorted(errors["M1_nf5", 6])
    print("\nMaximum relative error of the complete n_f = 5 coefficients (over angles and parameters)")
    print(f"{'beta':>7} {'M1 beta^4':>11} {'M1 beta^6':>11} {'M2 beta^4':>11} {'M2 beta^6':>11}")
    for b in betas:
        print(f"{float(b):7.3f}" + "".join(f" {errors[q, o][b]:11.2e}"
                                           for q in ("M1_nf5", "M2_nf5") for o in (4, 6)))
    print("\nFitted error slopes d ln(max error)/d ln(beta)")
    for q in ("M1_nf5", "M2_nf5"):
        for o in (4, 6):
            x = np.log([float(b) for b in betas]); y = np.log([errors[q, o][b] for b in betas])
            print(f"  {q} through beta^{o}: {np.polyfit(x, y, 1)[0]:.3f}")

    RESULTS.mkdir(exist_ok=True)
    slopes = {}
    for q in ("M1_nf5", "M2_nf5"):
        for o in (4, 6):
            x = np.log([float(b) for b in betas]); y = np.log([errors[q, o][b] for b in betas])
            slopes[f"{q}_order{o}"] = float(np.polyfit(x, y, 1)[0])
    json.dump({"records": count, "reproduction": float(worst_reproduction), "slopes": slopes,
               "max_error": {f"{q}_order{o}": {str(b): errors[q, o][b] for b in betas}
                             for (q, o) in errors}},
              open(RESULTS / "accuracy.json", "w"), indent=1)
    with open(RESULTS / "accuracy.csv", "w") as out:
        out.write("beta,M1_order4,M1_order6,M2_order4,M2_order6\n")
        for b in betas:
            out.write(f"{b}," + ",".join(f"{errors[q, o][b]:.6e}" for q in ("M1_nf5", "M2_nf5") for o in (4, 6)) + "\n")
    write_table(errors)
    print("\nwrote results/accuracy.csv and manuscript/tables/accuracy.tex")
    if worst_reproduction > REPRODUCTION_TOLERANCE:
        sys.exit(f"FAILED: the series evaluation differs from the stored values by "
                 f"{mp.nstr(worst_reproduction, 3)} (tolerance {REPRODUCTION_TOLERANCE:g})")
    print(f"Reproduction check passed (tolerance {REPRODUCTION_TOLERANCE:g}).")


SHOWN = [Fraction(1, 100), Fraction(1, 20), Fraction(1, 10), Fraction(1, 5), Fraction(3, 10)]


def write_table(errors):
    """Single-column table (Table III); the numbers come from the comparison above or results/accuracy.csv."""
    with open(TABLES / "accuracy.tex", "w") as out:
        out.write(TABLE_HEAD)
        for b in SHOWN:
            cells = [_latex(errors[q, o][b]) for q in ("M1_nf5", "M2_nf5") for o in (4, 6)]
            out.write(f"{float(b):.2f} & " + " & ".join(cells) + " \\\\\n")
        out.write("\\end{tabular}\n\\end{ruledtabular}\n\\end{table}\n")


def table_only():
    """Rewrite accuracy.tex from results/accuracy.csv without repeating the comparison."""
    errors = defaultdict(dict)
    for row in csv.DictReader(open(RESULTS / "accuracy.csv")):
        b = Fraction(row["beta"])
        for q, col in (("M1_nf5", "M1"), ("M2_nf5", "M2")):
            for o in (4, 6):
                errors[q, o][b] = float(row[f"{col}_order{o}"])
    write_table(errors)
    print("wrote manuscript/tables/accuracy.tex from results/accuracy.csv")


TABLE_HEAD = r"""\begin{table}[t]
\centering
\caption{Maximum relative errors of the complete Type II coefficients on the validation grid. At each value of $\beta_A$, the maximum is taken over the \AnglePointCount{} angles and \ParameterSetCount{} parameter sets of Table~\ref{tab:validation-grid}.}
\label{tab:AA-threshold-accuracy}
\small
\setlength{\tabcolsep}{3pt}
\begin{ruledtabular}
\begin{tabular}{ccccc}
 & \multicolumn{2}{c}{$\mathcal{M}_1$ through} & \multicolumn{2}{c}{$\mathcal{M}_2$ through} \\
$\beta_A$ & $\beta_A^4$ & $\beta_A^6$ & $\beta_A^4$ & $\beta_A^6$ \\
\hline
"""


def _latex(x):
    exponent = math.floor(math.log10(x))
    mantissa = x / 10**exponent
    return f"${mantissa:.2f}\\times 10^{{{exponent}}}$"


if __name__ == "__main__":
    table_only() if "--table-only" in sys.argv else main()
