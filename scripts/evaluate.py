#!/usr/bin/env python3
"""Evaluate M_1 and M_2 from the threshold series at one phase-space point.

Example:
    python3 scripts/evaluate.py --beta 1/10 --costheta 1/3

Output is the two-loop Type-II finite remainder with couplings and Wilson
coefficients stripped, in the conventions of arXiv:1909.08993.
"""

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import mpmath as mp

from ggaa.series import ThresholdSeries


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--beta", required=True, help="beta_A, e.g. 0.1 or 1/10")
    parser.add_argument("--costheta", required=True, help="cos(theta) in the partonic frame")
    parser.add_argument("--mA", default="125", help="pseudoscalar mass in GeV (default 125)")
    parser.add_argument("--muR", default="125", help="renormalization scale in GeV (default 125)")
    parser.add_argument("--Nc", default="3")
    parser.add_argument("--nf", default="5")
    parser.add_argument("--order", type=int, default=6, choices=(4, 6))
    parser.add_argument("--digits", type=int, default=30, help="working precision")
    args = parser.parse_args()

    mp.mp.dps = args.digits
    series = ThresholdSeries()
    for i in (1, 2):
        value = series(i, args.beta, args.costheta, mA=args.mA, muR=args.muR,
                       Nc=args.Nc, nf=args.nf, order=args.order)
        sign = "+" if value.imag >= 0 else "-"
        print(f"M{i} = {mp.nstr(value.real, args.digits - 5)} {sign} {mp.nstr(abs(value.imag), args.digits - 5)} i")


if __name__ == "__main__":
    main()
