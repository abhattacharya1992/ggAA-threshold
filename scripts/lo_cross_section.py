#!/usr/bin/env python3
"""Leading-order partonic cross section of the Type-II amplitude for g g -> A A.

In the effective theory of Eq. (3) of the paper, whose only coupling of A to gluons is
-(1/8) C_G Phi_A O_G, the two Type-II diagrams are the complete tree-level amplitude.
Prints sigma_hat and d sigma_hat / d cos(theta) and writes results/lo_cross_section.{csv,pdf}.

Example:
    python3 scripts/lo_cross_section.py --mA 125 --alpha-s 0.118 --cot-beta 1
"""

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp

from ggaa import benchmark
from ggaa.crosssection import dsigma_dcos, sigma
from ggaa.paths import RESULTS


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mA", type=float, default=float(benchmark.M_A), help="pseudoscalar mass in GeV")
    parser.add_argument("--alpha-s", type=float, default=benchmark.ALPHA_S)
    parser.add_argument("--cot-beta", type=float, default=float(benchmark.COT_BETA))
    parser.add_argument("--Nc", type=int, default=benchmark.N_C)
    args = parser.parse_args()
    mp.mp.dps = 30
    opts = dict(mA=args.mA, alpha_s=args.alpha_s, cot_beta=args.cot_beta, Nc=args.Nc)

    betas = [mp.mpf(k) / 100 for k in (1, 2, 5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 95)]
    RESULTS.mkdir(exist_ok=True)
    rows = []
    print(f"LO partonic cross section, m_A = {args.mA} GeV, alpha_s = {args.alpha_s}, cot(beta_H) = {args.cot_beta}")
    print(f"{'beta_A':>7} {'sqrt(s) [GeV]':>14} {'sigma_hat [pb]':>16}")
    for b in betas:
        rs = 2 * args.mA / mp.sqrt(1 - b**2)
        value = sigma(b, **opts)
        rows.append((b, rs, value))
        print(f"{float(b):7.2f} {float(rs):14.2f} {mp.nstr(value, 6):>16}")
    with open(RESULTS / "lo_cross_section.csv", "w") as out:
        out.write("beta,sqrt_s_GeV,sigma_hat_pb,mA_GeV,alpha_s,cot_beta,Nc\n")
        for b, rs, v in rows:
            out.write(f"{float(b)},{float(rs):.6f},{mp.nstr(v, 12)},{args.mA},{args.alpha_s},{args.cot_beta},{args.Nc}\n")

    fig, (left, right) = plt.subplots(1, 2, figsize=(7.15, 2.8))
    left.plot([float(r[1]) for r in rows], [float(r[2]) for r in rows], color="#8B2E3C", marker="o")
    left.set_xlabel(r"$\sqrt{\hat s}$ [GeV]")
    left.set_ylabel(r"$\hat\sigma_{\rm LO}$ [pb]")
    left.set_yscale("log")
    cos = [mp.mpf(k) / 20 for k in range(-19, 20)]
    for b, colour in ((mp.mpf("0.1"), "#202020"), (mp.mpf("0.5"), "#8B2E3C"), (mp.mpf("0.9"), "#1F7775")):
        right.plot([float(c) for c in cos], [float(dsigma_dcos(b, c, **opts)) for c in cos],
                   color=colour, label=rf"$\beta_A={float(b)}$")
    right.set_xlabel(r"$\cos\theta$")
    right.set_ylabel(r"$d\hat\sigma_{\rm LO}/d\cos\theta$ [pb]")
    right.set_yscale("log")
    right.legend(frameon=False)
    for ax in (left, right):
        ax.tick_params(direction="in", which="both", top=True, right=True)
    fig.tight_layout()
    fig.savefig(RESULTS / "lo_cross_section.pdf")
    print("\nwrote results/lo_cross_section.csv and results/lo_cross_section.pdf")


if __name__ == "__main__":
    main()
