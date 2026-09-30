#!/usr/bin/env python3
"""Independent evaluation of the exact two-loop remainders.

The published files Finite_2Loop_P1/P2 (run inputs/fetch_inputs.sh first) are evaluated
with ggaa.twoloop and ggaa.li22, which share no code with the Wolfram Language
implementation that produced the reference column of data/validation/validation_1080.csv.
The complete n_f = 5 coefficients are compared at points covering all three parameter
sets, and the flavour decomposition at the first point.  Each value is computed at two
working precisions to show that it is stable.

At cos(theta) = 0 (t = u) the published expressions contain removable 0/0 terms, so
they cannot be evaluated there numerically without taking a limit; the points below
avoid it.

Runtime: about 1.5 min (M_1) and 3 min (M_2) per evaluation, 45 min in total.
Use --points and --skip-flavour to run a subset.
"""

import argparse
import json
import pathlib
import sys
import time
from fractions import Fraction

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import mpmath as mp

from ggaa import twoloop
from ggaa.normalization import FILE_FACTOR
from ggaa.paths import INPUTS, RESULTS
from ggaa.validation import records, relative_error

POINTS = [  # (m_A, mu_R, N_c), beta, cos(theta)
    ((125, 125, 3), "1/10", "1/3"),
    ((125, 125, 3), "3/10", "-2/3"),
    ((125, 250, 3), "1/5", "-1/3"),
    ((150, 75, 4), "1/20", "2/3"),
]
PRECISIONS = (120, 160, 220, 300)
CACHE = RESULTS / "exact_amplitude_cache.json"   # finished evaluations, so an interrupted run can resume
LI22_STORE = RESULTS / "li22_store.json"          # Li22 values, for the same purpose


def exact(q):
    q = Fraction(q)
    return mp.mpf(q.numerator) / q.denominator


def evaluate(i, beta, costh, params, nf, digits):
    key = f"M{i} {params} beta={beta} cos={costh} nf={nf} digits={digits}"
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    if key in cache:
        re, im = cache[key]
        with mp.workdps(digits):
            return mp.mpc(re, im)
    mA, muR, Nc = params
    with mp.workdps(digits):
        value = twoloop.remainder(i, 2, exact(beta), exact(costh), mA, muR, Nc, nf, store=LI22_STORE) / FILE_FACTOR
        cache[key] = [mp.nstr(value.real, digits), mp.nstr(value.imag, digits)]
    RESULTS.mkdir(exist_ok=True)
    CACHE.write_text(json.dumps(cache, indent=1))
    return value


def stable_value(i, beta, costh, params, nf=5):
    """Raise the working precision until two successive values agree to 30 digits."""
    previous = None
    for digits in PRECISIONS:
        value = evaluate(i, beta, costh, params, nf, digits)
        if previous is not None and relative_error(previous, value) < mp.mpf(10) ** -30:
            return value, digits
        previous = value
    raise RuntimeError(f"M{i} at beta={beta}, cos={costh} not stable up to {PRECISIONS[-1]} digits")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--points", type=int, nargs="*", default=list(range(len(POINTS))),
                        help="indices into POINTS (default: all)")
    parser.add_argument("--form-factors", type=int, nargs="*", default=[1, 2])
    parser.add_argument("--skip-flavour", action="store_true")
    args = parser.parse_args()
    mp.mp.dps = 60   # for combining and comparing values; each evaluation sets its own precision
    if not (INPUTS / "Finite_2Loop_P1").exists():
        sys.exit("run inputs/fetch_inputs.sh first")
    table = {(r["quantity"], (r["mA"], r["muR"], r["Nc"]), r["beta"], r["costheta"]): r["exact"]
             for r in records()}
    worst = 0
    print("Complete n_f = 5 coefficients: independent evaluation against the validation table")
    for params, beta, costh in (POINTS[k] for k in args.points):
        for i in args.form_factors:
            start = time.time()
            value, digits = stable_value(i, beta, costh, params)
            reference = table[f"M{i}_nf5", params, Fraction(beta), Fraction(costh)]
            difference = relative_error(value, reference)
            worst = max(worst, difference)
            print(f"  M{i}  (mA, muR, Nc) = {params}  beta = {beta:>6}  cos = {costh:>4}:  "
                  f"stable at {digits} digits, relative difference {mp.nstr(difference, 3):>9}"
                  f"   ({time.time() - start:.0f} s)", flush=True)

    if args.skip_flavour:
        print(f"\nLargest relative difference: {mp.nstr(worst, 3)}")
        return
    params, beta, costh = POINTS[0]
    print(f"\nFlavour coefficients at (mA, muR, Nc) = {params}, beta = {beta}, cos = {costh} (160 digits)")
    for i in args.form_factors:
        v0, v1, v2 = (evaluate(i, beta, costh, params, nf, 160) for nf in (0, 1, 2))
        components = {0: v0, 2: (v2 - 2 * v1 + v0) / 2}
        components[1] = v1 - v0 - components[2]
        for k in (0, 1, 2):
            reference = table[f"M{i}_nf{k}", params, Fraction(beta), Fraction(costh)]
            difference = relative_error(components[k], reference)
            worst = max(worst, difference)
            print(f"  M{i}, n_f^{k}: relative difference {mp.nstr(difference, 3)}", flush=True)

    print(f"\nLargest relative difference: {mp.nstr(worst, 3)}")
    if worst > 1e-25:
        sys.exit("independent evaluation does not reproduce the validation table")


if __name__ == "__main__":
    main()
