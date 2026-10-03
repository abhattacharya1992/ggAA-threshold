#!/usr/bin/env python3
"""Figures 2 and 6 of the paper (files fig5_amplitudes, fig6_cross_section).

  fig5_amplitudes      |a_s^n M_i^(n)| / m_A^2 at tree level, one loop and two loops against beta_A
                       (the helicity amplitudes are C_G^2 M_1 for ++ and C_G^2 M_2 for +-)
  fig6_cross_section   LO partonic cross section and the finite virtual corrections a_s V1,
                       a_s^2 V2 (read from results/virtual_correction.csv, see `make virtual`)

Needs the published one-loop files (inputs/fetch_inputs.sh).
"""

import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp

from make_figures import BLACK, BLUE, GREY, RED, WIDE, finish_axis, headroom, save, style
from ggaa import benchmark
from ggaa.born import born
from ggaa.crosssection import sigma
from ggaa.oneloop import OneLoop
from ggaa.paths import INPUTS, RESULTS
from ggaa.series import ThresholdSeries

M_A = benchmark.M_A
COS = mp.mpf(benchmark.COS_THETA.numerator) / benchmark.COS_THETA.denominator


def fig_amplitudes():
    mp.mp.dps = 80
    one_loop, two_loop = OneLoop(), ThresholdSeries()
    betas = [mp.mpf(10) ** (e / mp.mpf(4)) for e in range(-16, -1)] + [mp.mpf("0.2"), mp.mpf("0.3")]
    fig, axes = plt.subplots(1, 2, figsize=(WIDE, 2.7))
    for i, ax in zip((1, 2), axes):
        a_s = benchmark.ALPHA_S / (4 * mp.pi)
        curves = {
            "tree": [abs(born(i, b, COS, M_A)) / M_A**2 for b in betas],
            r"$a_s\times$ one loop": [a_s * abs(one_loop(i, b, COS)) / M_A**2 for b in betas],
            r"$a_s^2\times$ two loops": [a_s**2 * abs(two_loop(i, b, COS)) / M_A**2 for b in betas],
        }
        x = [float(b) for b in betas]
        for (label, values), colour, marker, ls in zip(curves.items(), (GREY, BLACK, RED), ("o", "s", "^"), (":", "-", "--")):
            ax.loglog(x, [float(v) for v in values], color=colour, marker=marker, mfc=colour, mec=colour, ms=4.5, mew=0.8, ls=ls, lw=0.65, zorder=3, label=label)
        ax.set_xlabel(r"$\beta_A$")
        ax.set_ylabel(rf"$|a_s^n\,\mathcal{{M}}_{i}^{{(n)}}|/m_A^2$")
        headroom(ax, 0.25, log=True)
        finish_axis(ax, ("(a)" if i == 1 else "(b)") + rf"  $\mathcal{{M}}_{i}$" + (",  $J_z=0$" if i == 1 else ",  $|J_z|=2$"))
    axes[1].legend(loc="lower right")
    fig.tight_layout()
    save(fig, "fig5_amplitudes")


def fig_cross_section():
    mp.mp.dps = 30
    fig, (left, right) = plt.subplots(1, 2, figsize=(WIDE, 2.7))
    betas = [mp.mpf(10) ** (e / mp.mpf(4)) for e in range(-12, 0)] + [mp.mpf("0.9")]
    left.loglog([float(b) for b in betas], [float(sigma(b, mA=M_A, alpha_s=benchmark.ALPHA_S, cot_beta=benchmark.COT_BETA)) for b in betas],
                color=BLACK, marker="o", mfc=BLACK, mec=BLACK, ms=4.5, mew=0.8, lw=0.65, zorder=3)
    left.set_xlabel(r"$\beta_A$")
    left.set_ylabel(r"$\hat\sigma_{\rm LO}$ [pb]")
    headroom(left, 0.12, log=True)
    finish_axis(left, "(a)")
    with open(RESULTS / "virtual_correction.csv", newline="") as f:
    	rows = list(csv.DictReader(f))
    a_s = benchmark.ALPHA_S / (4 * float(mp.pi))
    x = [float(r["beta"]) for r in rows]
    right.semilogx( x, [a_s * float(r["V1"]) for r in rows], color=BLACK, marker="o", mfc=BLACK, mec=BLACK, ms=4.5, mew=0.8, lw=0.65, zorder=3, label=r"$a_s V_1$")
    right.semilogx( x, [a_s**2 * float(r["V2"]) for r in rows], color=RED, marker="s", mfc=RED, mec=RED, ms=4.5, mew=0.8, ls="--", lw=0.65, zorder=3, label=r"$a_s^2 V_2$")
    right.set_xlabel(r"$\beta_A$")
    right.set_ylabel(r"$a_s^n V_n$")
    headroom(right, 0.3)
    right.legend(loc="lower left")
    finish_axis(right, "(b)")
    fig.tight_layout()
    save(fig, "fig6_cross_section")


def main():
    style()
    if not (INPUTS / "Finite_1Loop_P1").exists():
        sys.exit("run inputs/fetch_inputs.sh first")
    fig_amplitudes()
    if (RESULTS / "virtual_correction.csv").exists():
        fig_cross_section()
    else:
        print("results/virtual_correction.csv missing: run `make virtual` for Figure 6")


if __name__ == "__main__":
    main()
