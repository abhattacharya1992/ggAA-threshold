#!/usr/bin/env python3
"""Manuscript tables generated from the data and results files.

  validation_grid.tex        data/validation/validation_1080.csv
  forward_angles.tex         data/validation/forward_angles.csv, results/accuracy.json
  mb_constants.tex           data/constants/mb_constants.csv, results/mb_constants.json
  leading_terms.tex          results/leading_logs.json
  timing.tex                 data/timing/timing_benchmark.csv
  computing_environment.tex  data/timing/timing_run.json

accuracy.tex and leading_logs.tex are written by compare_validation.py and
check_leading_logs.py.
"""

import csv
import json
import re
import pathlib
import sys
from collections import defaultdict
from fractions import Fraction

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import sympy as sp

from ggaa.paths import DATA, RESULTS

TABLES = DATA.parent / "manuscript" / "tables"


def sci(x):
    mantissa, exponent = f"{float(x):.2e}".split("e")
    return f"${mantissa}\\times10^{{{int(exponent)}}}$"


def forward_angles():
    rows = list(csv.DictReader(open(DATA / "validation" / "forward_angles.csv")))
    table = defaultdict(dict)
    for r in rows:
        if r["costheta"].startswith("-"):
            continue                       # equal to +99/100 by Bose symmetry (checked in the CSV)
        table[(r["quantity"][:2], r["costheta"])][r["beta"]] = (r["relative_error4"], r["relative_error6"])
    betas = sorted({r["beta"] for r in rows}, key=Fraction)
    positive = sorted({r["costheta"] for r in rows if not r["costheta"].startswith("-")}, key=Fraction)
    mirrored = sorted({r["costheta"] for r in rows if r["costheta"].startswith("-")}, key=Fraction)
    for c in mirrored:                      # Bose symmetry: the mirrored rows must repeat the positive ones
        for r in rows:
            if r["costheta"] == c:
                twin = next(t for t in rows if t["quantity"] == r["quantity"] and t["beta"] == r["beta"] and t["costheta"] == c[1:])
                assert twin["relative_error6"] == r["relative_error6"]
    lines = [r"\begin{table}[t]", r"\centering",
             r"\caption{Relative truncation errors $E_i^{(6)}$ of the complete $\nf=\BenchNf$ coefficients towards the beam, "
             r"at $\mA=\muR=\BenchMA\,\mathrm{GeV}$ and $\Nc=\BenchNc$. The last rows give the largest value on the grid of "
             r"Table~\ref{tab:validation-grid}." + (" The values at $c=" + ", ".join(f"{float(Fraction(c)):g}" for c in mirrored) +
             "$ equal those at $-c$." if mirrored else "") + "}",
             r"\label{tab:forward-angles}", r"\begin{ruledtabular}", r"\setlength{\tabcolsep}{3pt}", r"\begin{tabular}{ll" + "c" * len(betas) + "}",
             r" & $c$ & " + " & ".join(f"$\\bA={float(Fraction(b)):g}$" for b in betas) + r" \\", r"\hline"]
    for i in ("M1", "M2"):
        for c in positive:
            vals = " & ".join(sci(table[(i, c)][b][1]) for b in betas)
            lines.append(f"$\\mathcal M_{i[1]}$ & ${float(Fraction(c)):g}$ & {vals} \\\\")
    accuracy = json.load(open(RESULTS / "accuracy.json"))["max_error"]
    for i in ("M1", "M2"):
        grid = accuracy[f"{i}_nf5_order6"]
        lines.append(f"$\\mathcal M_{i[1]}$ & grid max. & " + " & ".join(sci(grid[str(Fraction(b))]) for b in betas) + " \\\\")
    lines += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table}"]
    (TABLES / "forward_angles.tex").write_text("\n".join(lines) + "\n")
    print("wrote manuscript/tables/forward_angles.tex")


def mb_constants():
    rows = list(csv.DictReader(open(DATA / "constants" / "mb_constants.csv")))
    lines = [r"\begin{table}[t]", r"\centering",
             r"\caption{Reductions of the endpoint constants of Eq.~\eqref{eq:mb-definition}, found with PSLQ at \PSLQDigits{} digits; "
             r"the largest residual is \PSLQResidual.}",
             r"\label{tab:mb-constants}", r"\begin{ruledtabular}", r"\small", r"\begin{tabular}{cl}",
             r"$[p,q,a,w]$ & reduction \\", r"\hline"]
    for r in rows:
        expr = sp.sympify(r["reduction"], locals={"polylog": sp.polylog, "zeta": sp.zeta})
        tex = sp.latex(expr).replace(r"\log{\left(2 \right)}", r"\ln 2").replace(r"\operatorname{Li}_{4}\left(\frac{1}{2}\right)", r"\Li_4(\tfrac12)")
        tex = tex.replace(r"\zeta\left(3\right)", r"\zeta_3")
        tex = re.sub(r"\\ln 2\^\{(\d)\}", r"\\ln^{\1}2", tex)
        lines.append(f"$[{r['p']},{r['q']},{r['a']},{r['w']}]$ & ${tex}$ \\\\")
    lines += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table}"]
    (TABLES / "mb_constants.tex").write_text("\n".join(lines) + "\n")
    print("wrote manuscript/tables/mb_constants.tex")


def frac_tex(x):
    f = Fraction(x)
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


def validation_grid():
    """Single-column Table II; long value lists are broken into lines of a few entries."""
    rows = list(csv.DictReader(open(DATA / "validation" / "validation_1080.csv")))
    betas = sorted({Fraction(r["beta"]) for r in rows})
    angles = sorted({Fraction(r["costheta"]) for r in rows})
    params = sorted({(int(r["mA_GeV"]), int(r["muR_GeV"]), int(r["Nc"])) for r in rows}, key=lambda t: (t[2], t[0], t[1]))

    def block(label, items, chunk):
        parts = [items[k:k + chunk] for k in range(0, len(items), chunk)]
        return [(label if n == 0 else "") + " & " + ", ".join(part) + ("," if n < len(parts) - 1 else "") + r" \\"
                for n, part in enumerate(parts)]

    lines = [r"\begin{table}[t]", r"\centering",
             r"\caption{Validation grid. At each point the six flavour coefficients and the two complete $\nf=\BenchNf$ coefficients are compared.}",
             r"\label{tab:validation-grid}", r"\begin{ruledtabular}", r"\begin{tabular}{ll}",
             r"quantity & values \\", r"\hline"]
    lines += block(r"$\bA$", [f"${frac_tex(b)}$" for b in betas], 3)
    lines += block(r"$c$", [f"${frac_tex(c)}$" for c in angles], 5)
    lines += block(r"$(\mA/\mathrm{GeV},\,\muR/\mathrm{GeV},\,\Nc)$", [f"$({m},{mu},{nc})$" for m, mu, nc in params], 1)
    lines += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table}"]
    (TABLES / "validation_grid.tex").write_text("\n".join(lines) + "\n")
    print("wrote manuscript/tables/validation_grid.tex")


def leading_terms():
    inventory = json.load(open(RESULTS / "leading_logs.json"))["inventory"]
    lines = [r"\begin{table}[t]", r"\centering",
             r"\caption{Leading terms of the threshold expansion of the two-loop Type II remainders. The second column lists the inverse powers with "
             r"nonzero coefficients in the construction files; they cancel identically after the reductions of "
             r"Sec.~\ref{sec:cancellations}, as do all odd powers. The last column gives the powers of $\mathcal L$ at leading power.}",
             r"\label{tab:leading-threshold-inventory}", r"\begin{ruledtabular}", r"\begin{tabular}{lccc}",
             r" & inverse powers & leading & powers \\", r"quantity & in construction & power & of $\mathcal L$ \\", r"\hline"]
    for i in (1, 2):
        for k in (0, 1, 2):
            item = inventory[f"M{i}_nf{k}"]
            inverse = ", ".join(f"$\\bA^{{{n}}}$" for n in item["inverse_powers_in_construction"]) or "--"
            lines.append(f"$\\Mflav{{{i}}}{{{k}}}$ & {inverse} & $\\bA^{{{item['leading_power']}}}$ & "
                         f"${','.join(str(p) for p in item['calL_powers'])}$" + r" \\")
    lines += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table}"]
    (TABLES / "leading_terms.tex").write_text("\n".join(lines) + "\n")
    print("wrote manuscript/tables/leading_terms.tex")


def timing():
    rows = list(csv.DictReader(open(DATA / "timing" / "timing_benchmark.csv")))
    lines = [r"\begin{table}[t]", r"\centering",
             r"\caption{Evaluation time for the two complete Type II coefficients at the benchmark parameters, "
             r"$c=\TimingCos$ and \TimingPrecision-digit working precision, for the unexpanded remainders (exact) "
             r"and the series. First calls start with empty caches; repeated times are medians of \TimingRepeats{} "
             r"evaluations at the same point. Starting the kernel and loading the files are excluded.}",
             r"\label{tab:timing}", r"\begin{ruledtabular}", r"\footnotesize", r"\setlength{\tabcolsep}{2.6pt}",
             r"\begin{tabular}{lcccccc}",
             r"& \multicolumn{3}{c}{first call [s]} & \multicolumn{3}{c}{repeated [s]} \\",
             r"$\bA$ & exact & series & ratio & exact & series & ratio \\", r"\hline"]
    for r in rows:
        assert r['CosTheta'] == rows[0]['CosTheta']
        lines.append(f"${r['Beta']}$ & {float(r['ExactFirstSeconds']):.2f} & {float(r['SeriesFirstSeconds']):.2f} & "
                     f"{float(r['FirstCallSpeedup']):.2f} & {float(r['ExactWarmMedianSeconds']):.4f} & "
                     f"{float(r['SeriesWarmMedianSeconds']):.4f} & {float(r['WarmSpeedup']):.2f}" + r" \\")
    lines += [r"\end{tabular}", r"\end{ruledtabular}", r"\end{table}"]
    (TABLES / "timing.tex").write_text("\n".join(lines) + "\n")
    print("wrote manuscript/tables/timing.tex")


def environment():
    run = json.load(open(DATA / "timing" / "timing_run.json"))
    env = run["environment"]
    text = ("The measurements used Wolfram Language " + env["wolfram"].strip('"').split(" (")[0] +
            f" under {env['os']} on an {env['cpu'].replace(' 16-Core Processor', ' processor')} with {env['ram_gib']}~GiB of memory." + "\n")
    (TABLES / "computing_environment.tex").write_text(text)
    print("wrote manuscript/tables/computing_environment.tex")


if __name__ == "__main__":
    validation_grid()
    forward_angles()
    mb_constants()
    leading_terms()
    timing()
    environment()
