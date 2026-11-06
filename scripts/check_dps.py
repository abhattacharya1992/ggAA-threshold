#!/usr/bin/env python3
"""The one-loop ln p_T term against the double-parton-scattering (DPS) singularity.

Gaunt and Stirling (JHEP 06 (2011) 048) factorize the DPS-singular part of a one-loop
amplitude into collinear splitting amplitudes of the incoming gluons and the hard
amplitudes gg -> A.  For g -> gg splitting with momentum fraction x, and hard amplitudes
of opposite sign for ++ and -- (pseudoscalar), the helicity sum gives

    F_A(x) = 1/(x(1-x)) - 2 x (1-x)        (a scalar would give + instead of -).

In the forward limit at fixed beta_A, x = (1+beta_A)/2.  The script compares the slope
d Re(M_1^(1),fin)/d ln(1-c), i.e. the coefficient of ln p_T^2, with pi N_c F_A(x) m_A^2.
At threshold F_A(1/2) = 7/2, which is the coefficient 7 pi N_c of calL = ln(2 p_T/m_A) + ...
Writes results/dps_check.txt.
"""

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import mpmath as mp

from ggaa.oneloop import OneLoop
from ggaa.paths import INPUTS, RESULTS


def main():
    if not (INPUTS / "Finite_1Loop_P1").exists():
        sys.exit("run inputs/fetch_inputs.sh first")
    mp.mp.dps = 300
    amp = OneLoop()
    m2 = 125**2
    F_A = lambda x: 1 / (x * (1 - x)) - 2 * x * (1 - x)
    F_S = lambda x: 1 / (x * (1 - x)) + 2 * x * (1 - x)
    lines = [f"{'N_c':>3} {'n_f':>3} {'beta':>6} {'slope (exact one loop)':>24} {'pi N_c F_A(x)':>20} {'rel. diff':>10} {'pi N_c F_S(x)':>14}"]
    worst = 0
    record = {"points": []}
    for Nc, nf in ((3, 5), (3, 0), (4, 2)):
        for b in ("0.001", "0.1", "0.3", "0.5", "0.7"):
            b = mp.mpf(b)
            v = [amp(1, b, 1 - mp.mpf(10) ** -e, Nc=Nc, nf=nf) / m2 for e in (14, 18)]
            slope = -(v[1].real - v[0].real) / (4 * mp.log(10))
            x = (1 + b) / 2
            pred = mp.pi * Nc * F_A(x)
            worst = max(worst, abs(slope / pred - 1))
            record["points"].append({"Nc": Nc, "nf": nf, "beta": float(b), "slope": float(slope),
                                     "F_A": float(F_A(x)), "F_S": float(F_S(x))})
            lines.append(f"{Nc:>3} {nf:>3} {float(b):6.3f} {mp.nstr(slope, 15):>24} {mp.nstr(pred, 15):>20} "
                         f"{mp.nstr(abs(slope / pred - 1), 2):>10} {mp.nstr(mp.pi * Nc * F_S(x), 8):>14}")
    lines.append(f"\nlargest relative difference to pi N_c F_A(x): {mp.nstr(worst, 2)}")
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "dps_check.txt").write_text("\n".join(lines) + "\n")
    record["max_relative_difference"] = float(worst)
    json.dump(record, open(RESULTS / "dps_check.json", "w"), indent=1)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
