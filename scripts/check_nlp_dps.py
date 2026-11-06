#!/usr/bin/env python3
"""DPS logarithms beyond leading power in the two-loop threshold series.

Checks, for general Nc and nf and through beta^6:
  (a) the coefficient of calL^3 in M1 at every order, divided by its beta^0 value;
  (b) the forward limit c -> 1 of the two-loop M1 logarithms, written in the
      Gram-determinant variable G = (1/2) ln[(mA^4 - t u)/mA^4] = calL - ln(1 - beta^2),
      against the DPS prediction C1 * [-(4/3) Nc G^3 - 2 beta0 G^2],
      C1 = 2 pi Nc F_A(x) mA^2, x = (1 + beta)/2;
  (c) the forward limit of the calL^2 coefficient of M2.
Run from the repository root:  python3 scripts/check_nlp_dps.py
"""
import pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))
import sympy as sp
from ggaa import wolfram as W
from ggaa.paths import SERIES

calL, G, nf, b = sp.symbols("calL G n_f beta")
ORDERS = (0, 2, 4, 6)

def series_in_calL(i):
    comps = [[W.reduce(e) for e in W.read(SERIES / f"M{i}_nf{k}.m")] for k in (0, 1, 2)]
    shift = W.calL_shift()
    return [sp.expand(sum(nf**k * comps[k][j] for k in (0, 1, 2)).subs(W.L, calL - shift))
            for j in range(4)]

def forward(e):
    e = sp.expand(e.subs(W.LMU, 0))
    return sp.expand(sp.limit(e, W.costh, 1)) if e.has(W.costh) else e

def main():
    M1, M2 = series_in_calL(1), series_in_calL(2)

    print("(a) M1, coefficient of calL^3 relative to beta^0 (all angles):")
    lead = M1[0].coeff(calL, 3)
    xh = (1 + b) / 2
    fa = sp.series(1 / (xh * (1 - xh)) - 2 * xh * (1 - xh), b, 0, 8).removeO()
    failures = []
    for j, n in enumerate(ORDERS):
        ratio = sp.factor(sp.cancel(M1[j].coeff(calL, 3) / lead))
        expected = fa.coeff(b, n) / fa.coeff(b, 0)
        print(f"    beta^{n}: {ratio}   (F_A expansion: {expected})")
        if sp.simplify(ratio - expected) != 0:
            failures.append(f"(a) beta^{n}")

    print("\n(b) M1, forward limit in G, minus the DPS prediction:")
    total = sum(b**n * M1[j] for j, n in enumerate(ORDERS))
    log1mb2 = sp.series(sp.log(1 - b**2), b, 0, 8).removeO()
    total = sp.expand(sp.series(sp.expand(total.subs(calL, G + log1mb2)), b, 0, 7).removeO())
    x = (1 + b) / 2
    C1 = sp.series(2 * sp.pi * W.nc * (1 / (x * (1 - x)) - 2 * x * (1 - x)), b, 0, 7).removeO()
    beta0 = (11 * W.nc - 2 * nf) / 3
    ok = True
    for n in ORDERS:
        for p, factor in ((3, -sp.Rational(4, 3) * W.nc), (2, -2 * beta0)):
            got = forward(total.coeff(b, n).coeff(G, p))
            pred = sp.expand(factor * C1.coeff(b, n) * W.m2 * W.PI / sp.pi)
            d = sp.simplify(got - pred)
            ok &= d == 0
            print(f"    beta^{n} G^{p}: {sp.factor(got)}   difference: {d}")
    print(f"    all match: {ok}")
    if not ok:
        failures.append("(b)")

    print("\n(c) M2, forward coefficient of calL^2 / (-4 Nc) per order, in units of pi Nc m2 / 2:")
    guess = sp.series(8 / (1 - b**2) - 7 - b**2, b, 0, 8).removeO()
    for j, n in enumerate(ORDERS):
        k1 = sp.simplify(forward(M2[j].coeff(calL, 2)) / (-4 * W.nc) / (W.PI * W.nc * W.m2 / 2))
        print(f"    beta^{n}: {k1}   (8/(1-beta^2) - 7 - beta^2: {guess.coeff(b, n)})")
        if sp.simplify(k1 - guess.coeff(b, n)) != 0:
            failures.append(f"(c) beta^{n}")
    if failures:
        raise SystemExit("FAILED: " + ", ".join(failures))
    print("\nAll checks passed.")

if __name__ == "__main__":
    main()
