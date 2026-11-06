#!/usr/bin/env python3
"""Logarithms of p_T in the finite virtual corrections V1, V2 at fixed angle.

With M_1^(0) -> 2 i m_A^2 and the leading-power forms of the one- and two-loop remainders,
the ln p_T terms of the amplitude are absorptive relative to the tree amplitude.  One expects
V1 -> (pi/2) A1 + Im A0 without logarithms, and V2 = a + b L + c L^2 with
L = ln(2 p_T/m_A) and c = -(7/4) pi^2 N_c^2 (from -14 pi^2 N_c^2 in 2 Re M0* M2 and
+(49/4) pi^2 N_c^2 in |M1|^2).  The script fits V2 at three very small velocities and writes results/rate_logs.json.

It repeats the fit at c = 0, 1/3, 2/3 (mu_R = m_A) and at c = 1/3 with mu_R = 2 m_A.  At leading
power v0 and v1 do not depend on the angle, and v1 = (7 pi N_c/2) Re A10 + 7 pi^2 N_c beta0 does not
depend on mu_R (Re A10 is scale independent); the scale dependence of V2 sits in v0.  With v1, v2 < 0,
V2 has its maximum at l* = -v1/(2 v2).
"""

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import mpmath as mp

from ggaa import benchmark
from ggaa.crosssection import VirtualCorrections
from ggaa.paths import INPUTS, RESULTS


def main():
    if not (INPUTS / "Finite_1Loop_P1").exists():
        sys.exit("run inputs/fetch_inputs.sh first")
    mp.mp.dps = 110
    V = VirtualCorrections()
    c = mp.mpf(benchmark.COS_THETA.numerator) / benchmark.COS_THETA.denominator
    Nc = benchmark.N_C
    rows = []
    for e in (6, 7, 8):
        b = mp.mpf(10) ** -e
        v1, v2, _ = V.at(b, c, Nc=Nc, nf=benchmark.N_F)
        L = mp.log(2 * b * mp.sqrt(1 - c**2) / mp.sqrt(1 - b**2))
        rows.append((e, L, v1, v2))
        print(f"beta = 1e-{e}: L = {mp.nstr(L, 8)}, V1 = {mp.nstr(v1, 12)}, V2 = {mp.nstr(v2, 12)}")
    A = mp.matrix([[1, r[1], r[1] ** 2] for r in rows])
    a, bcoef, ccoef = mp.lu_solve(A, mp.matrix([r[3] for r in rows]))
    expected = -mp.mpf(7) / 4 * mp.pi**2 * Nc**2
    print(f"fit: c = {mp.nstr(ccoef, 12)}, expected -(7/4) pi^2 N_c^2 = {mp.nstr(expected, 12)}")
    def fit(cc, muR):
        pts = []
        for e in (6, 7, 8):
            bb = mp.mpf(10) ** -e
            _, w2, _ = V.at(bb, cc, muR=muR, Nc=Nc, nf=benchmark.N_F)
            pts.append((mp.log(2 * bb * mp.sqrt(1 - cc**2) / mp.sqrt(1 - bb**2)), w2))
        return mp.lu_solve(mp.matrix([[1, q[0], q[0] ** 2] for q in pts]), mp.matrix([q[1] for q in pts]))

    angles = [mp.mpf(0), mp.mpf(1) / 3, mp.mpf(2) / 3]
    by_angle = [fit(cc, benchmark.MU_R) for cc in angles]
    at_2mu = fit(c, 2 * benchmark.MU_R)
    spread = max(max(abs(f[k] - by_angle[0][k]) / abs(by_angle[0][k]) for f in by_angle) for k in (0, 1))
    one = json.load(open(RESULTS / "one_loop.json"))
    re_a10 = mp.mpf([r for r in one["A0_scale"] if int(r["muR"]) == int(benchmark.MU_R)][0]["A0_re"]) \
        if isinstance(one["A0_scale"], list) else None
    beta0 = mp.mpf(11 * Nc - 2 * benchmark.N_F) / 3
    v1_formula = 7 * mp.pi * Nc / 2 * re_a10 + 7 * mp.pi**2 * Nc * beta0
    ell_peak = -bcoef / (2 * ccoef)
    beta_peak = mp.exp(ell_peak) / (2 * mp.sqrt(1 - c**2))
    print(f"v0, v1 at c = 0, 1/3, 2/3: " + "; ".join(f"{mp.nstr(f[0], 10)}, {mp.nstr(f[1], 10)}" for f in by_angle))
    print(f"largest relative spread of v0, v1 over the angles: {mp.nstr(spread, 3)}")
    print(f"c = 1/3, mu_R = 2 m_A: v0 = {mp.nstr(at_2mu[0], 10)}, v1 = {mp.nstr(at_2mu[1], 10)}")
    print(f"v1 formula (7 pi Nc/2) Re A10 + 7 pi^2 Nc beta0 = {mp.nstr(v1_formula, 10)} (fit {mp.nstr(bcoef, 10)})")
    print(f"maximum of V2 at l* = {mp.nstr(ell_peak, 6)}, beta* = {mp.nstr(beta_peak, 4)} for c = 1/3")
    RESULTS.mkdir(exist_ok=True)
    json.dump({"cos": float(c), "Nc": Nc, "betas": [f"1e-{r[0]}" for r in rows],
               "V1": [float(r[2]) for r in rows], "fit_a": float(a), "fit_b": float(bcoef), "fit_c": float(ccoef),
               "expected_c": float(expected), "relative_difference": float(abs(ccoef / expected - 1)),
               "angle_fits": [{"cos": float(cc), "v0": float(f[0]), "v1": float(f[1]), "v2": float(f[2])} for cc, f in zip(angles, by_angle)],
               "angle_spread_v0_v1": float(spread),
               "scale_fit_2muR": {"cos": float(c), "v0": float(at_2mu[0]), "v1": float(at_2mu[1]), "v2": float(at_2mu[2])},
               "v1_scale_rel_diff": float(abs(at_2mu[1] / bcoef - 1)),
               "v1_formula": float(v1_formula), "v1_formula_rel_diff": float(abs(v1_formula / bcoef - 1)),
               "ell_peak": float(ell_peak), "beta_peak": float(beta_peak)},
              open(RESULTS / "rate_logs.json", "w"), indent=1)


if __name__ == "__main__":
    main()
