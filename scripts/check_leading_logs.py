#!/usr/bin/env python3
"""Leading power of the threshold series in terms of one logarithm.

With calL = ln(beta) + [ln(1 - c) + ln(1 + c)]/2 + ln 2 + i pi/2 = (1/2) ln[-4 beta^2 (1 - c^2)]
the script checks, for every order and flavour coefficient, that no ln(1 - c),
ln(1 + c) or Li2((1 + c)/2) remains once ln(beta) is traded for calL, and prints
the beta^0 coefficients of the complete M_1 and M_2 for general n_f.
"""

import pathlib
import json
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import sympy as sp

from ggaa import benchmark
from ggaa import wolfram as W
from ggaa.paths import CONSTRUCTION, RESULTS, SERIES, TABLES

calL, nf = sp.symbols("calL n_f")


def load(i, k):
    return [W.reduce(c) for c in W.read(SERIES / f"M{i}_nf{k}.m")]


def in_calL(expr):
    return sp.expand(expr.subs(W.L, calL - W.calL_shift()))


def numeric(expr, Nc=benchmark.N_C, n_f=benchmark.N_F):
    """Value of a beta^0 coefficient divided by m_A^2 at mu_R = m_A."""
    values = {W.m2: 1, W.mu2: 1, W.LMU: 0, W.LN2: sp.log(2), W.PI: sp.pi, W.Z3: sp.zeta(3),
              W.LI4: sp.polylog(4, sp.Rational(1, 2)), W.nc: Nc, nf: n_f}
    return complex(sp.N(expr.subs(values), 20))


# exact beta^0 coefficients of calL^k (units of pi m_A^2, mu_R = m_A), checked against the reduction
EXACT_FORMS = {(1, 3): -sp.Rational(28, 3) * W.nc**2, (1, 2): -sp.Rational(14, 3) * W.nc * (11 * W.nc - 2 * nf),
               (2, 2): -2 * W.nc**2, (2, 1): -sp.Rational(4, 3) * W.nc * (4 * W.nc - nf)}
EXACT = {(1, 3): r"$-\tfrac{28}{3}\pi N_c^2$", (1, 2): r"$-\tfrac{14}{3}\pi N_c(11N_c-2n_f)$",
         (1, 1): r"Eq.~\eqref{eq:M1-L1}", (2, 2): r"$-2\pi N_c^2$", (2, 1): r"$-\tfrac{4}{3}\pi N_c(4N_c-n_f)$"}


def write_table(table):
    """manuscript/tables/leading_logs.tex: the beta^0 coefficients of calL^k."""
    rows = []
    for i in (1, 2):
        for power in sorted(table[i], reverse=True):
            value = numeric(table[i][power])
            exact_tex = EXACT.get((i, power), "--")
            if (i, power) in EXACT_FORMS:
                difference = sp.expand(table[i][power].subs(W.LMU, 0) / (W.PI * W.m2) - EXACT_FORMS[i, power])
                assert difference == 0, (i, power, difference)
            rows.append(f"$\\mathcal M_{i}$ & ${power}$ & {exact_tex} & ${_tex(value.real)}$ & ${_tex(value.imag)}$ \\\\")
    body = "\n".join(rows)
    text = TABLE.replace("BODY", body)
    target = TABLES / "leading_logs.tex"
    target.write_text(text)
    print(f"wrote {target}")


def _tex(x):
    if abs(x) < 1e-12:
        return "0"
    mantissa, exponent = f"{x:.5e}".split("e")
    return f"{mantissa}\\times10^{{{int(exponent)}}}"


TABLE = r"""\begin{table*}[!htb]
\centering
\caption{Coefficients of the leading power term, $\mathcal M_i|_{\beta_A^0}=m_A^2\sum_k \widehat A_{ik}\,\mathcal L^k$, with $\mathcal L$ of Eq.~\eqref{eq:calL}. They do not depend on the scattering angle. Exact forms for $\mu_R=m_A$ and general $N_c$, $n_f$; numerical values for $N_c=\BenchNc$, $n_f=\BenchNf$, $\mu_R=m_A$.}
\label{tab:leading-log-benchmark}
\begin{ruledtabular}
\begin{tabular}{cclrr}
 & $k$ & exact & $\operatorname{Re}\widehat A_{ik}$ & $\operatorname{Im}\widehat A_{ik}$ \\
\hline
BODY
\end{tabular}
\end{ruledtabular}
\end{table*}
"""


def main():
    clean = True
    table = {}
    inventory = {}
    checked = 0
    for i in (1, 2):
        components = [load(i, k) for k in (0, 1, 2)]
        for n, index in zip((0, 2, 4, 6), range(4)):
            for k in (0, 1, 2):
                terms = W.collect(in_calL(components[k][index]).subs(calL, W.L))
                left = [m for m in terms if m.has(W.LM) or m.has(W.LP) or m.has(W.LI2)]
                angular = [m for m in terms if terms[m].has(W.costh)]
                expr = in_calL(components[k][index])
                even = sp.simplify(expr - expr.subs(W.costh, -W.costh)) == 0
                checked += 1
                clean &= not left and even
                print(f"M{i} nf^{k} beta^{n}: {len(terms):3d} terms, angular functions outside calL: "
                      f"{len(left)}, rational angular dependence: {'yes' if angular else 'no'}, "
                      f"even in cos(theta): {even}")
        complete = sum(nf**k * components[k][0] for k in (0, 1, 2))
        leading = sp.collect(in_calL(complete), calL, evaluate=False)
        print(f"\nM{i} at beta^0, in units of pi m_A^2, for general n_f:")
        degree = max(sp.degree(m, calL) for m in leading)
        for power in range(degree, 0, -1):
            coefficient = sp.factor(leading[calL**power] / (W.PI * W.m2))
            print(f"   calL^{power}: {coefficient}")
        print()
        table[i] = {power: leading.get(calL**power if power else sp.S.One, 0) for power in range(degree + 1)}
        for k in (0, 1, 2):
            inverse = [n for n in (-4, -3, -2, -1)
                       if (CONSTRUCTION / f"M{i}_nf{k}_beta_{'m' + str(-n)}.wl").exists()
                       and W.read(CONSTRUCTION / f"M{i}_nf{k}_beta_{'m' + str(-n)}.wl") != 0]
            first = next(n for n, c in zip((0, 2, 4, 6), components[k]) if sp.expand(c) != 0)
            lead = sp.expand(in_calL(components[k][(0, 2, 4, 6).index(first)]))
            powers = sorted({int(sp.degree(m, calL)) for m in sp.Add.make_args(lead)})
            inventory[f"M{i}_nf{k}"] = {"inverse_powers_in_construction": inverse, "leading_power": first,
                                        "calL_powers": powers}
    write_table(table)
    record = {"angle_only_through_calL_and_even": bool(clean), "coefficients_checked": checked,
              "inventory": inventory,
              "exact_units_pi_mA2": {f"{i},{power}": sp.srepr(sp.factor(table[i][power] / (W.PI * W.m2)))
                                     for i in table for power in table[i] if power >= 1},
              "numeric_Nc3_nf5_muR_eq_mA": {f"{i},{power}": [numeric(table[i][power]).real, numeric(table[i][power]).imag]
                                           for i in table for power in table[i]}}
    json.dump(record, open(RESULTS / "leading_logs.json", "w"), indent=1)
    print("Angle enters only through calL, and every coefficient is even in cos(theta):", clean)
    if not clean:
        sys.exit(1)


if __name__ == "__main__":
    main()
