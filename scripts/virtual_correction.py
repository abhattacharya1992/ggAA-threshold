#!/usr/bin/env python3
"""Finite virtual corrections to the Type-II amplitude, relative to LO.

    sum|M|^2 = C_G^4 sum_i |M_i^(0)|^2 [ 1 + a_s V1 + a_s^2 V2 ],   a_s = alpha_s / (4 pi),

weighted with the LO angular distribution.  V1 uses the exact one-loop remainder; V2 uses
the two-loop threshold series and is therefore limited to beta_A <~ 0.3.  Both are finite
remainders after Catani's infrared subtraction and are scheme dependent.

All orders are in one normalization (ggaa.normalization).  V2 does not include the O_G-O_J
term proportional to the physical C_J, which depends on m_t.

Writes results/virtual_correction.{csv,pdf}.
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
from ggaa.crosssection import VirtualCorrections
from ggaa.paths import INPUTS, RESULTS


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alpha-s", type=float, default=benchmark.ALPHA_S)
    parser.add_argument("--nodes", type=int, default=64, help="Gauss-Legendre nodes in cos(theta)")
    parser.add_argument("--reference-nodes", type=int, default=128, help="nodes of the convergence check")
    args = parser.parse_args()
    if not (INPUTS / "Finite_1Loop_P1").exists():
        sys.exit("run inputs/fetch_inputs.sh first")
    mp.mp.dps = 60
    a_s = mp.mpf(args.alpha_s) / (4 * mp.pi)
    corrections = VirtualCorrections()

    betas = [mp.mpf(v) for v in ("0.001", "0.003", "0.01", "0.03", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3")]
    rows = []
    print(f"a_s = alpha_s/(4 pi) = {mp.nstr(a_s, 6)}, angle-averaged over {args.nodes} nodes")
    print(f"{'beta':>6} {'V1':>10} {'V2':>10} {'a_s V1':>8} {'a_s^2 V2':>9} {'1 + a_s V1 + a_s^2 V2':>22}")
    for b in betas:
        v1, v2 = corrections.angle_averaged(b, nodes=args.nodes)
        rows.append((b, v1, v2))
        print(f"{float(b):6.3f} {mp.nstr(v1, 7):>10} {mp.nstr(v2, 7):>10} {mp.nstr(a_s * v1, 4):>8} "
              f"{mp.nstr(a_s**2 * v2, 4):>9} {mp.nstr(1 + a_s * v1 + a_s**2 * v2, 4):>22}", flush=True)

    RESULTS.mkdir(exist_ok=True)
    with open(RESULTS / "virtual_correction.csv", "w") as out:
        out.write("beta,V1,V2\n")
        for b, v1, v2 in rows:
            out.write(f"{float(b)},{mp.nstr(v1, 12)},{mp.nstr(v2, 12)}\n")

    # convergence of the angular quadrature, quoted in Sec. V B
    check = ("0.01", "0.1", "0.3")
    dev1 = dev2 = 0.0
    for b in check:
        a1, a2 = corrections.angle_averaged(mp.mpf(b), nodes=args.nodes)
        r1, r2 = corrections.angle_averaged(mp.mpf(b), nodes=args.reference_nodes)
        dev1 = max(dev1, float(abs(a1 / r1 - 1))); dev2 = max(dev2, float(abs(a2 / r2 - 1)))
    print(f"quadrature: {args.nodes} against {args.reference_nodes} nodes at beta = {', '.join(check)}: "
          f"V1 changes by at most {dev1:.1e}, V2 by at most {dev2:.1e}")
    import json
    json.dump(dict(nodes=args.nodes, reference_nodes=args.reference_nodes, betas=list(check),
                   max_dev_V1=dev1, max_dev_V2=dev2), open(RESULTS / "virtual_quadrature.json", "w"), indent=1)

    fig, ax = plt.subplots(figsize=(3.6, 2.8))
    x = [float(r[0]) for r in rows]
    ax.semilogx(x, [float(a_s * r[1]) for r in rows], color="#202020", marker="o", label=r"$a_s V_1$")
    ax.semilogx(x, [float(a_s**2 * r[2]) for r in rows], color="#8B2E3C", marker="s", ls="--", label=r"$a_s^2 V_2$")
    ax.axhline(0, color="#6F6F6F", lw=0.7)
    ax.set_xlabel(r"$\beta_A$")
    ax.set_ylabel("correction relative to LO")
    ax.legend(frameon=False, fontsize=7)
    ax.tick_params(direction="in", which="both", top=True, right=True)
    fig.tight_layout()
    fig.savefig(RESULTS / "virtual_correction.pdf")
    print("\nwrote results/virtual_correction.csv and results/virtual_correction.pdf")


if __name__ == "__main__":
    main()
