#!/usr/bin/env python3
"""Size of the logarithmic terms of V2 at fixed angle.

Uses the leading-power fit V2 = v0 + v1 l + v2 l^2 of scripts/check_rate_logs.py
(results/rate_logs.json, c = 1/3) and the complete V2 from ggaa.crosssection at the same angle.
Writes results/log_share.json.
"""

import json
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import mpmath as mp

from ggaa.crosssection import VirtualCorrections
from ggaa.paths import RESULTS

BETAS = (0.01, 0.1)


def main():
    mp.mp.dps = 40
    fit = json.load(open(RESULTS / "rate_logs.json"))
    c = fit["cos"]
    v0, v1, v2 = fit["fit_a"], fit["fit_b"], fit["fit_c"]
    V = VirtualCorrections()
    out = {"cos": c, "v0": v0, "v1": v1, "v2": v2, "points": []}
    for b in BETAS:
        l = math.log(2 * b * math.sqrt(1 - c * c) / math.sqrt(1 - b * b))
        _, V2, _ = V.at(mp.mpf(b), mp.mpf(1) / 3)
        V2 = float(V2)
        row = {"beta": b, "l": l, "V2": V2, "leading_power": v0 + v1 * l + v2 * l * l,
               "log_share": (v1 * l + v2 * l * l) / V2, "square_share": v2 * l * l / V2}
        out["points"].append(row)
        print(f"beta = {b}: l = {l:.3f}, V2 = {V2:.1f} (leading power {row['leading_power']:.1f}); "
              f"(v1 l + v2 l^2)/V2 = {row['log_share']:.3f}, v2 l^2/V2 = {row['square_share']:.3f}")
    json.dump(out, open(RESULTS / "log_share.json", "w"), indent=1)


if __name__ == "__main__":
    main()
