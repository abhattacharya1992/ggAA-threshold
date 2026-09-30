#!/usr/bin/env python3
"""Velocity below which the two-loop term of M_1 exceeds the one-loop term (Fig. 4, Sec. IV D).

Towards threshold |a_s^2 M_1^(2)| grows as L^3 and |a_s M_1^(1)| as L, so the two-loop term overtakes
the one-loop term at small beta_A.  The crossing is located by bisection in ln(beta_A) at the
benchmark angle, mu_R = m_A and alpha_s of the benchmark, with the one-loop remainder evaluated
exactly (including the scale term) and the two-loop coefficient from the series through beta_A^6.
Writes results/amplitude_crossover.json.
"""

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import mpmath as mp

from ggaa import benchmark
from ggaa.oneloop import OneLoop
from ggaa.paths import RESULTS
from ggaa.series import ThresholdSeries


def main():
    mp.mp.dps = 80
    one, two = OneLoop(), ThresholdSeries()
    c = mp.mpf(benchmark.COS_THETA.numerator) / benchmark.COS_THETA.denominator
    a_s = benchmark.ALPHA_S / (4 * mp.pi)

    def excess(lnb):
        b = mp.exp(lnb)
        return mp.log(abs(a_s**2 * two(1, b, c)) / abs(a_s * one(1, b, c)))

    lo, hi = mp.log(mp.mpf("1e-4")), mp.log(mp.mpf("0.3"))
    assert excess(lo) > 0 > excess(hi), "no single crossing in the plotted range"
    for _ in range(60):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if excess(mid) > 0 else (lo, mid)
    b = mp.exp((lo + hi) / 2)
    print(f"|a_s^2 M1^(2)| = |a_s M1^(1)| at beta_A = {mp.nstr(b, 6)} (c = {mp.nstr(c, 4)}, mu_R = m_A)")
    json.dump({"cos": float(c), "alpha_s": benchmark.ALPHA_S, "beta_crossover": float(b)},
              open(RESULTS / "amplitude_crossover.json", "w"), indent=1)


if __name__ == "__main__":
    main()
