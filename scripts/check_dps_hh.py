#!/usr/bin/env python3
"""One-loop gg -> HH through two CP-even effective vertices: the forward logarithm (Sec. 4.3, Fig. 4).

Reads data/dps_hh/dps_hh.csv and data/dps_hh/hh_dps.log (see data/dps_hh/README.md) and checks
- that the HH forward coefficient of ln pT^2, in units of pi N_c m^2 and multiplied by the sign of the
  HH tree amplitude at threshold relative to AA (-1), equals F_S((1+beta)/2), recomputed here;
- that the AA coefficient from the same run equals F_A;
- that the |J_z| = 2 amplitude has no logarithm;
- the reproduction of the 2019 AA result and the double-pole check recorded in the log;
- the threshold ratio of the two coefficients, each in units of its tree amplitude, against
  F_S/F_A at the velocity of the threshold point (beta = 1e-6), 9/7 at beta = 0.
Writes results/dps_hh.json.
"""

import csv
import json
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

from ggaa.paths import RESULTS

DATA = pathlib.Path(__file__).resolve().parents[1] / "data" / "dps_hh"


def FA(x):
    return 1 / (x * (1 - x)) - 2 * x * (1 - x)


def FS(x):
    return 1 / (x * (1 - x)) + 2 * x * (1 - x)


def main():
    rows = list(csv.DictReader(open(DATA / "dps_hh.csv")))
    points = []
    for r in rows:
        b = float(r["beta"]); x = (1 + b) / 2
        aa, hh = float(r["AA"]), float(r["HH"])
        points.append({"beta": b, "AA": aa, "HH": hh, "HH_tree_sign": -hh,
                       "dev_AA": aa / FA(x) - 1, "dev_HH": -hh / FS(x) - 1,
                       "M2_over_M1": float(r["HH_M2/HH_M1"])})
    log = (DATA / "hh_dps.log").read_text()
    # [4] reproduction of the archived 2019 AA result: "branch +1: 0``40." = zero to 40 digits
    m = re.search(r"branch \+1: 0``(\d+)", log)
    reproduce_digits = int(m.group(1)) if m else None
    # [5] double poles against Catani's operator: every ratio printed as 1.`15.
    poles = re.findall(r"\(1/eps\^2 pole\)/\(-8 Nc tree\) = ([0-9.]+)`(\d+)", log)
    pole_ok = len(poles) == 4 and all(abs(float(v) - 1) < 1e-12 for v, _ in poles)
    pole_digits = min(int(d) for _, d in poles) if poles else None
    tree = re.search(r"HH tree at threshold, M1\(0\)/mA\^2: (-?[0-9.]+)", log)
    tree_hh = float(tree.group(1))            # coefficient of i
    # [6] threshold slopes at beta = 1e-6, in units of pi N_c
    thr = re.search(r"threshold slope/\(pi Nc\):\s+AA ([-0-9.]+)`\S*\s+HH ([-0-9.]+)`", log)
    s_aa, s_hh = float(thr.group(1)), float(thr.group(2))
    b0 = 1e-6; x0 = (1 + b0) / 2
    tree_aa = 2.0 / (1 - b0**2)               # AA tree i s/2 at the same velocity (the log's HH tree is at beta = 1e-6 too)
    ratio = (s_hh / tree_hh) / (s_aa / tree_aa)
    ratio_dev = ratio / (FS(x0) / FA(x0)) - 1
    out = {
        "points": points,
        "beta_min": min(p["beta"] for p in points), "beta_max": max(p["beta"] for p in points),
        "max_dev_HH": max(abs(p["dev_HH"]) for p in points),
        "max_dev_AA": max(abs(p["dev_AA"]) for p in points),
        "max_M2_over_M1": max(p["M2_over_M1"] for p in points),
        "aa_reproduction_digits": reproduce_digits,
        "double_poles_ok": pole_ok, "double_pole_digits": pole_digits,
        "tree_HH_threshold_over_i": tree_hh, "tree_AA_threshold_over_i": tree_aa,
        "threshold_ratio_in_tree_units": ratio, "threshold_ratio_dev_from_FS_over_FA": ratio_dev,
    }
    problems = []
    if out["max_dev_HH"] > 1e-12: problems.append("HH does not follow F_S")
    if out["max_dev_AA"] > 1e-12: problems.append("AA does not follow F_A")
    if out["max_M2_over_M1"] > 1e-10: problems.append("|J_z| = 2 amplitude has a logarithm")
    if not reproduce_digits or reproduce_digits < 30: problems.append("2019 AA result not reproduced")
    if not pole_ok: problems.append("double-pole check missing or failed")
    if abs(tree_hh + tree_aa) > 1e-9: problems.append("unexpected HH tree at threshold")
    if abs(ratio_dev) > 1e-12: problems.append("threshold ratio differs from F_S/F_A")
    print(f"HH forward coefficient vs F_S: largest deviation {out['max_dev_HH']:.1e} "
          f"for {out['beta_min']} <= beta <= {out['beta_max']} (AA vs F_A in the same run: {out['max_dev_AA']:.1e})")
    print(f"|M2 slope|/|M1 slope| for HH: at most {out['max_M2_over_M1']:.1e}")
    print(f"2019 AA result reproduced to {reproduce_digits} digits; double poles as Catani's operator "
          f"({pole_digits} digits); HH tree at threshold {tree_hh:+.0f}i m^2")
    print(f"threshold ratio HH/AA in tree units: {ratio:.15f} (F_S/F_A at beta = 1e-6: "
          f"{FS(x0) / FA(x0):.15f}; 9/7 = {9 / 7:.15f}); deviation {ratio_dev:.1e}")
    RESULTS.mkdir(exist_ok=True)
    json.dump(out, open(RESULTS / "dps_hh.json", "w"), indent=1)
    if problems:
        sys.exit("PROBLEMS: " + "; ".join(problems))


if __name__ == "__main__":
    main()
