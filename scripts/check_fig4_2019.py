#!/usr/bin/env python3
"""Which normalization does Fig. 4 of arXiv:1909.08993 use?

Fig. 4 there plots M^(2),fin_i / M^(0)_i against x for several cos(theta).  This script
evaluates that ratio with the two-loop remainder in the package normalization (the
published files divided by 64) and shows it on two axes: left, the consistent ratio;
right, the same curve times 64, which is what the figure shows if the published two-loop
file was divided by the unscaled tree amplitude.  Uses the sixth-order series, so it
covers x >= 0.54 (beta_A <= 0.3).  Writes results/fig4_2019_check.{txt,pdf}.
"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp

from make_figures import INK, TEAL, WINE, finish_axis, style
from ggaa.born import born
from ggaa.normalization import FILE_FACTOR
from ggaa.paths import RESULTS
from ggaa.series import ThresholdSeries

ANGLES = (mp.mpf("0.2"), mp.mpf("0.5"), mp.mpf("0.8"))


def main():
    mp.mp.dps = 40
    series = ThresholdSeries()
    xs = [mp.mpf(k) / 100 for k in range(54, 100, 2)]
    beta = lambda x: (1 - x) / (1 + x)
    rows = []
    for c in ANGLES:
        for x in xs:
            r = [series(i, beta(x), c) / born(i, beta(x), c) for i in (1, 2)]
            rows.append((c, x, r))
    lines = ["ratio M^(2),fin_i / M^(0)_i (m_A = mu_R = 125 GeV, N_c = 3, n_f = 5)",
             "consistent normalization; multiply by 64 for the published-file normalization",
             f"{'cos':>5} {'x':>5} {'Re M1':>12} {'Im M1':>12} {'Re M2':>12} {'Im M2':>12}"]
    for c, x, r in rows:
        if int(x * 100) % 10 == 4:
            lines.append(f"{float(c):5.2f} {float(x):5.2f} " + " ".join(f"{float(v):12.1f}" for z in r for v in (z.real, z.imag)))
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "fig4_2019_check.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    style()
    fig, axes = plt.subplots(1, 2, figsize=(7.15, 2.9))
    for i, ax in zip((0, 1), axes):
        for c, colour in zip(ANGLES, (INK, WINE, TEAL)):
            sel = [(float(x), r[i]) for cc, x, r in rows if cc == c]
            ax.plot([p[0] for p in sel], [float(p[1].real) for p in sel], color=colour, ls="-", label=rf"Re, $c={float(c)}$")
            ax.plot([p[0] for p in sel], [float(p[1].imag) for p in sel], color=colour, ls="--", label=rf"Im, $c={float(c)}$")
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(rf"$\mathcal{{M}}_{i+1}^{{(2),\rm fin}}/\mathcal{{M}}_{i+1}^{{(0)}}$ (consistent)")
        ax.secondary_yaxis("right", functions=(lambda y: y * FILE_FACTOR, lambda y: y / FILE_FACTOR)).set_ylabel(r"$\times 64$ (published files)")
        finish_axis(ax, rf"$\mathcal{{M}}_{i+1}$")
    axes[0].legend(frameon=False, fontsize=6)
    fig.tight_layout()
    fig.savefig(RESULTS / "fig4_2019_check.pdf")
    print("wrote results/fig4_2019_check.txt and results/fig4_2019_check.pdf")


if __name__ == "__main__":
    main()
