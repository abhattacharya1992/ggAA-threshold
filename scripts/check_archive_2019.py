#!/usr/bin/env python3
"""Compare the 2019 archive outputs (data/provenance_2019/archive) with the published files.

1. One loop.  The archive remainders Finite1L1/2 are ratios to the tree amplitude.  They depend
   on mu_R as 2 beta0 ln mu_R^2 and agree with the corrected remainders of ggaa.oneloop
   (published + 2 beta0 ln(mu_R^2/GeV^2) M^(0)); the published files equal them at mu_R = 1 GeV.
2. Two loops.  The archive run (no factor 64, O_J term switched off with CJ = 0) is evaluated at
   its phase-space point x = 1/11, z = 3/11 (beta = 5/6, cos = -29/30).  The published P1
   remainder is compared with 64 x (archive x tree) + (8/6) Z_GJ^(1) [mu^-eps G(eps)]_0, the
   O_G-O_J counterterm of the published pipeline, with G the one-loop O_G-O_J amplitude ggjAA11.
   The published value needs 120 digits and takes about 3 minutes; it is read from
   results/archive_2019_published_P1.json unless --compute is given.
Writes results/archive_2019.json.
"""

import json
import lzma
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "python"))

import mpmath as mp

from ggaa.born import born
from ggaa.oneloop import OneLoop, _to_python
from ggaa.paths import RESULTS
from ggaa.twoloop import Evaluator, remainder

ARCHIVE = ROOT / "data" / "provenance_2019" / "archive"
CACHE = RESULTS / "archive_2019_published_P1.json"


def kinematics(b, c, mA):
    x = (1 - b) / (1 + b)
    y = (1 + b**2 - 2 * b * c) / (1 - b**2)
    z = (1 + b**2 + 2 * b * c) / (1 - b**2)
    M2 = mp.mpf(mA) ** 2
    s, t = M2 * (1 + x) ** 2 / x, -M2 * y
    return x, y, z, M2, s, t


def one_loop_check(record):
    one = OneLoop()
    beta0 = mp.mpf(23) / 3
    worst_fixed, shifts, offsets = 0, [], []
    for i in (1, 2):
        code = compile(_to_python(lzma.open(ARCHIVE / "one_loop" / f"Finite1L{i}.xz", "rt").read()), f"F1L{i}", "eval")
        for b, c in ((mp.mpf("0.3"), mp.mpf("0.4")), (mp.mpf("0.1"), mp.mpf(1) / 3)):
            x, y, z, M2, s, t = kinematics(b, c, 125)
            u = 2 * M2 - s - t

            def archive(muR):
                return eval(code, dict(mpf=mp.mpf, log=mp.log, pi=mp.pi, j=mp.mpc(0, 1),
                                       polylog=lambda n, a: mp.polylog(int(n), a), x=x, z=z, Y=1 + y, Z=1 + z,
                                       M2=M2, muR2=mp.mpf(muR) ** 2, N=mp.mpf(3), Nf=mp.mpf(5), PT=(u * t - M2**2) / s))

            M0 = born(i, b, c, 125)
            a125, a250 = archive(125), archive(250)
            fixed = one(i, b, c, mA=125, muR=125) / M0
            raw = one(i, b, c, mA=125, muR=125, scale_fix=False) / M0
            worst_fixed = max(worst_fixed, abs(a125 - fixed) / abs(fixed))
            shifts.append(abs((a250 - a125) / (2 * beta0 * mp.log(4)) - 1))
            offsets.append(abs((a125 - raw) / (2 * beta0 * mp.log(mp.mpf(125) ** 2)) - 1))
            print(f"   M{i} beta={mp.nstr(b, 2)} cos={mp.nstr(c, 3)}: archive - corrected = {mp.nstr(abs(a125 - fixed), 3)}, "
                  f"archive - published = {mp.nstr((a125 - raw).real, 12)} x M0")
    record["one_loop"] = {"max_rel_dev_from_corrected": mp.nstr(worst_fixed, 3),
                          "max_rel_dev_mu_shift": mp.nstr(max(shifts), 3),
                          "max_rel_dev_offset_2b0_ln_mu2": mp.nstr(max(offsets), 3)}


def two_loop_check(record, compute):
    mp.mp.dps = 40
    b, c = mp.mpf(5) / 6, mp.mpf(-29) / 30
    x, y, z, M2, s, t = kinematics(b, c, 125)
    assert abs(x - mp.mpf(1) / 11) < 1e-30 and abs(z - mp.mpf(3) / 11) < 1e-30
    text = lzma.open(ARCHIVE / "two_loop" / "Finite2L1.xz", "rt").read().replace("\\\n", "")
    sym = {"Y": 1 + y, "Z": 1 + z, "M2": M2, "muR2": M2, "N": mp.mpf(3), "Nf": mp.mpf(5), "Pi": +mp.pi,
           "I": mp.mpc(0, 1), "CJ": mp.mpf(0), "s": s, "t": t, "zeta2": mp.pi**2 / 6, "zeta3": mp.zeta(3)}
    arch = Evaluator(sym).evaluate(text)
    code = compile(_to_python((ARCHIVE / "two_loop" / "ggjAA11").read_text()), "gj", "eval")

    def G(e):
        return eval(code, dict(mpf=mp.mpf, log=mp.log, pi=mp.pi, j=mp.mpc(0, 1), zeta=lambda n: mp.zeta(int(n)),
                               x=x, z=z, ma2=M2, nf=mp.mpf(5), CJ=mp.mpf(1), eps=mp.mpf(e)))

    h = mp.mpf("1e-12")
    G0, G1 = G(0), (G(h) - G(-h)) / (2 * h)
    M0, CF = born(1, b, c, 125), mp.mpf(4) / 3
    oj = mp.mpf(8) / 6 * (-24 * CF) * (G1 - mp.log(M2) / 2 * G0)
    if compute or not CACHE.exists():
        mp.mp.dps = 120
        pub = remainder(1, 2, b, c, mA=125, muR=125)
        json.dump([mp.nstr(pub.real, 40), mp.nstr(pub.imag, 40)], open(CACHE, "w"))
        mp.mp.dps = 40
    pub = mp.mpc(*map(mp.mpf, json.load(open(CACHE))))
    pred = 64 * arch * M0 + oj
    re_dev = abs(pub.real - pred.real) / abs(pub.real)
    residual = (pub - pred) / M0
    print(f"   published P1 = {mp.nstr(pub, 22)}\n   64 x archive x tree + O_J = {mp.nstr(pred, 22)}")
    print(f"   real parts: relative difference {mp.nstr(re_dev, 3)}; remaining term / M0 = {mp.nstr(residual, 12)}"
          f" = {mp.nstr(residual.real / (32 * CF * G0 / M0).real, 12)} x 32 C_F G0/M0")
    record["two_loop"] = {"point": {"beta": "5/6", "cos": "-29/30", "x": "1/11", "z": "3/11"},
                          "archive_ratio": [mp.nstr(arch.real, 25), mp.nstr(arch.imag, 25)],
                          "published_P1": [mp.nstr(pub.real, 25), mp.nstr(pub.imag, 25)],
                          "real_part_rel_dev": mp.nstr(re_dev, 3),
                          "G0_over_M0": mp.nstr((G0 / M0).real, 20),
                          "residual_over_M0": mp.nstr(residual.real, 20),
                          "residual_over_32CF_G0": mp.nstr(residual.real / (32 * CF * G0 / M0).real, 20)}


def main():
    mp.mp.dps = 60
    record = {}
    print("1. One loop: archive remainders (ratio to tree) against the published and corrected files")
    one_loop_check(record)
    print("\n2. Two loops at x = 1/11, z = 3/11 (P1)")
    two_loop_check(record, "--compute" in sys.argv)
    json.dump(record, open(RESULTS / "archive_2019.json", "w"), indent=1)
    print("\nwrote results/archive_2019.json")


if __name__ == "__main__":
    main()
