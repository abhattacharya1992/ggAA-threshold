#!/usr/bin/env python3
"""Top-mass dependence of the Type II tree amplitude near threshold.

In the full theory the Type II tree amplitude consists of two top-quark triangles, each
coupling one A to an on-shell and a virtual gluon, connected by t- or u-channel gluon
exchange.  For AA its top-mass dependence is 2 [I_2(tau_A, lambda)]^2 per channel, with
tau_A = 4 m_t^2/m_A^2 and lambda = 4 m_t^2/t (or u), in terms of the A -> Z gamma loop
function I_2 (Baglio et al., arXiv:2303.05409, Eq. (30)); I_2 -> 1/2 for m_t -> infinity.
The script evaluates [I_2/I_2^HTL]^2 at threshold (t = u = -m_A^2).

For the DPS logarithm the effective vertices enter through the amplitude g g -> A with both
gluons on shell (the 2 -> 1 amplitudes of the Gaunt-Stirling factorisation, evaluated at zero
transverse momentum).  That is the limit lambda -> infinity of I_2, I_2 -> (tau/2) f(tau), so
the coefficient of the logarithm is multiplied by [tau f(tau)]^2 relative to the heavy-top limit.
Writes results/top_mass.json.
"""

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import mpmath as mp

from ggaa import benchmark
from ggaa.paths import RESULTS

MASSES = (125, 150, 170)


def f(tau):
    tau = mp.mpf(tau)
    if tau >= 1:
        return mp.asin(1 / mp.sqrt(tau)) ** 2
    r = mp.sqrt(1 - tau)
    return -(mp.log((1 + r) / (1 - r)) - mp.mpc(0, 1) * mp.pi) ** 2 / 4


def I2(tau, lam):
    return -tau * lam / (2 * (tau - lam)) * (f(tau) - f(lam))


def main():
    mp.mp.dps = 30
    mt = mp.mpf(benchmark.M_TOP)
    out = {"m_t": float(mt), "threshold_ratio": {}}
    for mA in MASSES:
        tau = 4 * mt**2 / mA**2
        ratio = (I2(tau, 4 * mt**2 / (-mA**2)) / mp.mpf(0.5)) ** 2
        assert abs(ratio.imag) < 1e-20
        out["threshold_ratio"][str(mA)] = float(ratio.real)
        dps = (tau * f(tau)) ** 2
        assert abs(mp.im(dps)) < 1e-20
        out.setdefault("dps_vertex_ratio", {})[str(mA)] = float(mp.re(dps))
        print(f"m_A = {mA} GeV: full/HTL at threshold = {mp.nstr(ratio.real, 6)}, DPS vertices [tau f(tau)]^2 = {mp.nstr(mp.re(dps), 6)}")
    json.dump(out, open(RESULTS / "top_mass.json", "w"), indent=1)


if __name__ == "__main__":
    main()
