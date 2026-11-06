#!/usr/bin/env python3
"""Check the two inputs of the construction of the n_f coefficient of M1.

src/construction/inputs/P1_nf1.wl is the coefficient of n_f in the published two-loop P1
remainder, and P1_nf1_compact.wl a compact form of it.  The construction reads both and
checks their SHA-256.  This script checks their content:

1. P1_nf1.wl, evaluated with ggaa.twoloop, against the coefficient of n_f of M1 in
   data/validation/validation_1080.csv, which was obtained from the published file
   independently of P1_nf1.wl, at one point of each parameter set (published normalisation,
   40 stored digits);
2. P1_nf1_compact.wl against P1_nf1.wl at two further points.

Writes results/nf1_inputs.txt.  Takes about half a minute.
"""

import csv
import hashlib
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "python"))

import mpmath as mp

from ggaa.paths import RESULTS
from ggaa.twoloop import Evaluator, symbols

INPUTS = ROOT / "src" / "construction" / "inputs"
SHA256 = {"P1_nf1.wl": "0da89f4e7a22e433478b3e810b13cd720cf6c5eb29e168f771d1acd200be3815",
          "P1_nf1_compact.wl": "af362f3a8757a4335570347ccbbd781f87d060df87d163ab123bf37ac21958f4"}
TABLE_POINTS = [("125", "125", "3", "1/10", "1/3"), ("125", "250", "3", "1/5", "-2/3"),
                ("150", "75", "4", "3/10", "2/3")]
COMPACT_POINTS = [("1/7", "-1/5", 125, 125, 3), ("2/9", "3/7", 150, 75, 4)]


def fraction(text):
    p, q = (int(v) for v in text.split("/"))
    return mp.mpf(p) / q


def value(text, beta, cos, mA, muR, Nc):
    return Evaluator(symbols(fraction(beta), fraction(cos), mA, muR, Nc, 0)).evaluate(text)


def main():
    mp.mp.dps = 60
    lines, ok = [], True
    for name, digest in SHA256.items():
        same = hashlib.sha256((INPUTS / name).read_bytes()).hexdigest() == digest
        ok &= same
        lines.append(f"{name}: SHA-256 {'as required by the construction' if same else 'DIFFERS'}")
    full = (INPUTS / "P1_nf1.wl").read_text().replace("\\\n", "")
    compact = (INPUTS / "P1_nf1_compact.wl").read_text().replace("\\\n", "").strip()
    table = {(r["quantity"], r["mA_GeV"], r["muR_GeV"], r["Nc"], r["beta"], r["costheta"]): r
             for r in csv.DictReader(open(ROOT / "data" / "validation" / "validation_1080.csv"))}
    for mA, muR, Nc, beta, cos in TABLE_POINTS:
        row = table[("M1_nf1", mA, muR, Nc, beta, cos)]
        ref = mp.mpc(mp.mpf(row["re_exact"]), mp.mpf(row["im_exact"]))
        diff = abs(value(full, beta, cos, int(mA), int(muR), int(Nc)) - ref) / abs(ref)
        ok &= diff < mp.mpf("1e-38")
        lines.append(f"P1_nf1.wl vs validation table, mA={mA} muR={muR} Nc={Nc} beta={beta} c={cos}: "
                     f"relative difference {mp.nstr(diff, 3)}")
    for beta, cos, mA, muR, Nc in COMPACT_POINTS:
        a, b = value(full, beta, cos, mA, muR, Nc), value(compact, beta, cos, mA, muR, Nc)
        diff = abs(a - b) / abs(a)
        ok &= diff < mp.mpf("1e-40")
        lines.append(f"P1_nf1_compact.wl vs P1_nf1.wl, mA={mA} muR={muR} Nc={Nc} beta={beta} c={cos}: "
                     f"relative difference {mp.nstr(diff, 3)}")
    lines.append("All checks passed." if ok else "FAILED")
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "nf1_inputs.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
