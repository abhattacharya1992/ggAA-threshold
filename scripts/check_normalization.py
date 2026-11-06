#!/usr/bin/env python3
"""Relative normalization of the tree, one-loop and two-loop amplitudes.

1. Born against BornSqr.log (n = 4): (|M_1^(0)|^2 + |M_2^(0)|^2) / BornSqr must be constant.
2. Born against the two Feynman diagrams evaluated with explicit polarization vectors:
   |A(++)| = |M_1^(0)| and |A(+-)| = |M_2^(0)|.  This fixes the absolute normalization
   (amplitude = C_G^2 M_i delta^{ab}).
3. Born against the one-loop remainder: the n_f term of M^(1),fin divided by M^(0) is the
   same for both tensor structures, and its imaginary part is -2 pi/3, the value produced by
   the infrared subtraction 2 I^(1) M^(0) with beta_0 = 11/3 N_c - 2/3 n_f.
4. Two-loop against one-loop: renormalization-group invariance of C_G^2 [M^(0) + a_s M^(1) +
   a_s^2 M^(2)] with C_G ~ a_s gives dM^(2)/d ln mu^2 = a beta_0 M^(1) + b M^(0) with a small
   integer a, 3 for renormalized amplitudes.  With the two-loop series in the package
   normalization (the published files divided by 64, ggaa.normalization) the fit gives a = 3;
   the published files themselves give 192.
"""

import json
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import mpmath as mp
import sympy as sp

from ggaa.born import born, helicity_amplitudes
from ggaa.oneloop import OneLoop
from ggaa.normalization import FILE_FACTOR
from ggaa.paths import DATA, INPUTS, RESULTS
from ggaa.series import ThresholdSeries

mp.mp.dps = 50


def born_sqr():
    text = (DATA / "born" / "original" / "BornSqr.log").read_text()
    body = text[text.index("BornSqrCF * (") + len("BornSqrCF * ("):text.rindex(");")]
    s, t, u, ma2, n = sp.symbols("s t u ma2 n")
    expr = sp.sympify(body.replace("^", "**").replace("\n", " "), locals=dict(s=s, t=t, u=u, ma2=ma2, n=n))
    return sp.lambdify((s, t, u, ma2), expr.subs(n, 4), "mpmath")


def main():
    points = [(mp.mpf("0.1"), mp.mpf(1) / 3), (mp.mpf("0.3"), mp.mpf(-2) / 3), (mp.mpf("0.6"), mp.mpf("0.2"))]
    m = 125

    record = {"file_factor": FILE_FACTOR, "bornsqr_ratio": [], "feynman_deviation": [], "nf_over_born": [], "rg_a": []}
    print("1. Born against BornSqr.log (n = 4)")
    bsq = born_sqr()
    for b, c in points:
        s = 4 * m**2 / (1 - b**2)
        t, u = m**2 - s / 2 * (1 - b * c), m**2 - s / 2 * (1 + b * c)
        ratio = (abs(born(1, b, c)) ** 2 + abs(born(2, b, c)) ** 2) / bsq(s, t, u, m**2)
        print(f"   beta={mp.nstr(b, 2)} cos={mp.nstr(c, 3):>6}: ratio = {mp.nstr(ratio, 15)}")
        record["bornsqr_ratio"].append(mp.nstr(ratio, 30))

    print("\n2. Born against the Feynman diagrams")
    for b, c in points:
        h = helicity_amplitudes(b, c, m)
        record["feynman_deviation"].append(max(abs(abs(h[1, 1]) / abs(complex(born(1, b, c))) - 1),
                                               abs(abs(h[1, -1]) / abs(complex(born(2, b, c))) - 1)))
        print(f"   beta={mp.nstr(b, 2)} cos={mp.nstr(c, 3):>6}: |A(++)|/|M1| = {abs(h[1, 1]) / abs(complex(born(1, b, c))):.12f}"
              f"   |A(+-)|/|M2| = {abs(h[1, -1]) / abs(complex(born(2, b, c))):.12f}")

    if not (INPUTS / "Finite_1Loop_P1").exists():
        sys.exit("\nrun inputs/fetch_inputs.sh for checks 3 and 4")
    one_loop = OneLoop()
    series = ThresholdSeries()
    print("\n3. n_f term of the one-loop remainder divided by the Born (expect equal values, Im = -2 pi/3 = -2.0944)")
    for b, c in ((mp.mpf(10) ** -6, mp.mpf(1) / 3), (mp.mpf("0.2"), mp.mpf(-2) / 3)):
        for i in (1, 2):
            with mp.workdps(80):
                r = (one_loop(i, b, c, nf=1) - one_loop(i, b, c, nf=0)) / born(i, b, c)
            print(f"   beta={mp.nstr(b, 2)} cos={mp.nstr(c, 3):>6} M{i}: {mp.nstr(r, 12)}")
            record["nf_over_born"].append({"beta": mp.nstr(b, 3), "cos": mp.nstr(c, 6), "i": i,
                                           "re": mp.nstr(r.real, 20), "im": mp.nstr(r.imag, 20)})

    print("\n4. Dimensional analysis: M/m_A^2 at fixed beta, c and mu_R/m_A must not change when m_A and mu_R are scaled together")
    record["scaling"] = []
    for b, c in ((mp.mpf("0.1"), mp.mpf(1) / 3), (mp.mpf("0.3"), mp.mpf(-2) / 3)):
        for i in (1, 2):
            with mp.workdps(60):
                raw = [one_loop(i, b, c, mA=mm, muR=mm, scale_fix=False) / mm**2 for mm in (125, 250)]
                fixed = [one_loop(i, b, c, mA=mm, muR=mm) / mm**2 for mm in (125, 250)]
                beta0 = mp.mpf(23) / 3
                expect = -2 * beta0 * mp.log(4) * born(i, b, c, 125) / 125**2
                dev_raw = abs((raw[1] - raw[0]) - expect) / abs(expect)
                dev_fixed = abs(fixed[1] - fixed[0]) / abs(fixed[0])
            print(f"   beta={mp.nstr(b, 2)} cos={mp.nstr(c, 3):>6} M{i}: published files change by -2 beta0 ln4 M0 "
                  f"(relative deviation {mp.nstr(dev_raw, 2)}); with the scale term the change is {mp.nstr(dev_fixed, 2)}")
            record["scaling"].append({"beta": mp.nstr(b, 3), "cos": mp.nstr(c, 6), "i": i,
                                      "raw_dev": mp.nstr(dev_raw, 5), "fixed_dev": mp.nstr(dev_fixed, 5)})

    print("\n5. Scale dependence of the two-loop series: dM2/dln(mu^2) = a beta0 M1 + b M0")
    print("   expect a = 3 and b = 2 beta1 - 8 H2g - 12 C_F n_f (Catani operators and O_G-O_J counterterm of the 2019 code)")
    series = ThresholdSeries()
    z2, z3 = mp.pi**2 / 6, mp.zeta(3)
    record["rg"] = []
    for (b, c) in ((mp.mpf("0.05"), mp.mpf(1) / 3), (mp.mpf("0.1"), mp.mpf(-2) / 3)):
        for Nc, nf, mm in ((3, 5, 125), (4, 2, 150), (3, 0, 300)):
            Ca, Cf = mp.mpf(Nc), mp.mpf(Nc**2 - 1) / (2 * Nc)
            beta0 = (11 * Ca - 2 * nf) / 3
            beta1 = mp.mpf(34) / 3 * Ca**2 - 2 * nf * Cf - mp.mpf(10) / 3 * Ca * nf
            H2g = (Ca**2 * (-mp.mpf(5) / 12 - mp.mpf(11) / 24 * z2 - z3 / 2) + Ca * nf * (mp.mpf(29) / 27 + z2 / 12)
                   - Cf * nf / 2 - nf**2 * mp.mpf(5) / 27) / 2
            predicted_b = 2 * beta1 - 8 * H2g - 12 * Cf * nf
            rows, rhs = [], []
            for i in (1, 2):
                v = [series(i, b, c, mA=mm, muR=mm * mp.e ** (l / 2), Nc=Nc, nf=nf) for l in (-1, 1)]
                rhs.append((v[1] - v[0]) / 2)
                rows.append([beta0 * one_loop(i, b, c, mA=mm, muR=mm, Nc=Nc, nf=nf), born(i, b, c, mm)])
            a, bcoef = mp.lu_solve(mp.matrix(rows), mp.matrix(rhs))
            v = [series(1, b, c, mA=mm, muR=mm * mp.e ** (l / 2), Nc=Nc, nf=nf) for l in (-2, 0, 2)]
            second = (v[2] - 2 * v[1] + v[0]) / 4 / (6 * beta0**2 * born(1, b, c, mm))
            print(f"   beta={mp.nstr(b, 2)} cos={mp.nstr(c, 3):>6} N_c={Nc} n_f={nf} m_A={mm}: a = {mp.nstr(a.real, 10)}, "
                  f"b = {mp.nstr(bcoef.real, 10)} (expected {mp.nstr(predicted_b, 10)}), d2M2/dL2 / (6 beta0^2 M0) = {mp.nstr(second.real, 10)}")
            record["rg"].append({"beta": mp.nstr(b, 3), "cos": mp.nstr(c, 6), "Nc": Nc, "nf": nf, "mA": mm,
                                 "a": mp.nstr(a.real, 20), "b": mp.nstr(bcoef.real, 20), "b_expected": mp.nstr(predicted_b, 20),
                                 "second": mp.nstr(second.real, 20)})
    json.dump(record, open(RESULTS / "normalization.json", "w"), indent=1)


if __name__ == "__main__":
    main()
