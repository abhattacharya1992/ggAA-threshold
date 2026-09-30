#!/usr/bin/env python3
"""One-loop Type-II remainders near threshold and in the forward direction.

Uses the published one-loop expressions (run inputs/fetch_inputs.sh first).
  1. Reproduces the reference values in data/validation/one_loop_reference.csv.
  2. Leading power: M = m_A^2 (A0 + A1 calL) + O(beta^2) with angle-independent
     A0, A1; prints A1 and M_2(beta -> 0) for several N_c, n_f and mu_R.
  3. Fixed beta, cos(theta) -> 1: M_1 grows like ln(1 - c), M_2 tends to a constant.
"""

import csv
import json
import pathlib
import sys
from fractions import Fraction

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import mpmath as mp

from ggaa.oneloop import OneLoop
from ggaa.paths import RESULTS, VALIDATION

M2 = 125**2


def exact(q):
    q = Fraction(q)
    return mp.mpf(q.numerator) / q.denominator


def calL_shift(c):
    return (mp.log(1 - c) + mp.log(1 + c)) / 2 + mp.log(2) + mp.mpc(0, 1) * mp.pi / 2


def main():
    amp = OneLoop()

    mp.mp.dps = 60
    print("1. Reference values")
    with open(VALIDATION / "one_loop_reference.csv") as handle:
        for row in csv.DictReader(handle):
            value = amp(int(row["form_factor"][1]), exact(row["beta"]), exact(row["costheta"]), scale_fix=False)
            reference = mp.mpc(row["re"], row["im"])
            print(f"   {row['form_factor']} beta={row['beta']:>5} cos={row['costheta']:>4}: "
                  f"relative difference {mp.nstr(abs(value - reference) / abs(reference), 3)}")

    mp.mp.dps = 400
    tiny = [mp.mpf(10) ** -40, mp.mpf(10) ** -50]
    record = {"fit_betas": [mp.nstr(b, 3) for b in tiny], "digits": mp.mp.dps, "stored_digits": 30, "angles": [], "scan": []}
    print("\n2. Leading power, M = m_A^2 (A0 + A1 calL)")
    for c in (mp.mpf(0), mp.mpf(1) / 3, mp.mpf(2) / 3, mp.mpf(-9) / 10):
        v1 = [amp(1, b, c) / M2 for b in tiny]
        slope = (v1[0] - v1[1]) / (mp.log(tiny[0]) - mp.log(tiny[1]))
        a0 = v1[0] - slope * (mp.log(tiny[0]) + calL_shift(c))
        m2_limit = amp(2, tiny[0], c) / M2
        print(f"   cos={mp.nstr(c, 3):>6}:  M1: A1 = {mp.nstr(slope.real, 12)}, A0 = {mp.nstr(a0, 10)};"
              f"  M2(beta->0)/m_A^2 = {mp.nstr(m2_limit.real, 12)}")
        record["angles"].append({"cos": mp.nstr(c, 6), "A1": mp.nstr(slope.real, 30), "A0_re": mp.nstr(a0.real, 30),
                                 "A0_im": mp.nstr(a0.imag, 30), "M2_limit": mp.nstr(m2_limit.real, 30)})
    print("   Colour, flavour and scale dependence (in units of pi m_A^2, cos = 1/3):")
    c = mp.mpf(1) / 3
    for Nc, nf, muR in ((3, 5, 125), (4, 5, 125), (3, 0, 125), (5, 2, 125), (3, 5, 250)):
        v1 = [amp(1, b, c, muR=muR, Nc=Nc, nf=nf) / M2 for b in tiny]
        slope = (v1[0] - v1[1]) / (mp.log(tiny[0]) - mp.log(tiny[1]))
        m2_limit = amp(2, tiny[0], c, muR=muR, Nc=Nc, nf=nf) / M2
        print(f"     N_c={Nc} n_f={nf} mu_R={muR:>3}:  A1 = {mp.nstr(slope.real / mp.pi, 10)} pi,"
              f"  M2(beta->0) = {mp.nstr(m2_limit.real / mp.pi, 10)} pi")
        record["scan"].append({"Nc": Nc, "nf": nf, "muR": muR, "A1_over_pi": mp.nstr(slope.real / mp.pi, 30),
                               "M2_limit_over_pi": mp.nstr(m2_limit.real / mp.pi, 30)})
    print("   Constant A0 against the scale (with the scale term of ggaa.oneloop; expect dA0/dln mu^2 = 4 i beta0):")
    record["A0_scale"] = []
    for muR in (125, 250):
        v1 = [amp(1, b, c, muR=muR) / M2 for b in tiny]
        slope = (v1[0] - v1[1]) / (mp.log(tiny[0]) - mp.log(tiny[1]))
        a0 = v1[0] - slope * (mp.log(tiny[0]) + calL_shift(c))
        print(f"     mu_R={muR}: A0 = {mp.nstr(a0, 12)}")
        record["A0_scale"].append({"muR": muR, "A0_re": mp.nstr(a0.real, 30), "A0_im": mp.nstr(a0.imag, 30)})
    json.dump(record, open(RESULTS / "one_loop.json", "w"), indent=1)

    mp.mp.dps = 300
    print("\n3. Fixed beta, cos(theta) -> 1 (values / m_A^2)")
    for beta in ("0.3", "0.1", "0.01"):
        for e in (2, 4, 6, 8, 12):
            c = 1 - mp.mpf(10) ** -e
            m1, m2 = amp(1, mp.mpf(beta), c) / M2, amp(2, mp.mpf(beta), c) / M2
            print(f"   beta={beta:>4}  1-c=1e-{e:<2}  M1 = {mp.nstr(m1, 10):>34}   M2 = {mp.nstr(m2, 10)}")


if __name__ == "__main__":
    main()
