#!/usr/bin/env python3
"""Figures 1, 3, 4 and 5 of the paper (files fig1-fig4), all from data/validation/validation_1080.csv.

  fig1_validation      truncation errors of the complete n_f = 5 coefficients
  fig2_flavour         |a_r| / |M_i(5)| for the terms n_f^0, 5 n_f^1, 25 n_f^2
  fig3_interference    magnitude deficit D_i and signed interference I_i
  fig4_ratio           modulus and phase of M_2 / M_1

Figures 2-4 use the unexpanded remainders at m_A = mu_R = 125 GeV, N_c = 3.
"""

import math
import pathlib
import sys
from collections import defaultdict
from fractions import Fraction

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ggaa import benchmark
from ggaa.paths import RESULTS, FIGURES
from ggaa.validation import records, relative_error

# Conventional journal style: black and the standard matplotlib colours, always paired with distinct
# line styles or open symbols (APS: curves must be distinguishable in grey scale as well).
BLACK, BLUE, RED, GREEN, GREY = "k", "#1f77b4", "#d62728", "#2ca02c", "0.45"
INK, WINE, TEAL, LIGHT = BLACK, RED, BLUE, "0.6"          # names used by the other figure scripts
BENCHMARK = (benchmark.M_A, benchmark.MU_R, benchmark.N_C)
BENCHMARK_ANGLE = benchmark.COS_THETA
WIDE = 7.0                                                  # two columns (17.8 cm); one column is 3.375 in


def style():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["STIXGeneral", "Times New Roman", "DejaVu Serif"],
        "mathtext.fontset": "stix", "font.size": 8, "axes.labelsize": 8.5,
        "legend.fontsize": 7.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "lines.linewidth": 1.0, "lines.markersize": 3.6, "lines.markeredgewidth": 0.8,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
        "xtick.major.size": 3.5, "ytick.major.size": 3.5, "xtick.minor.size": 2.0, "ytick.minor.size": 2.0,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "legend.frameon": False, "legend.handlelength": 2.2, "legend.borderaxespad": 0.6,
        "axes.grid": False, "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    })


def finish_axis(ax, label, loc=(0.04, 0.93)):
    """Boxed axes with inward major and minor ticks on all sides; panel label inside the frame."""
    from matplotlib.ticker import AutoMinorLocator
    for axis, scale in ((ax.xaxis, ax.get_xscale()), (ax.yaxis, ax.get_yscale())):
        if scale == "linear":
            axis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="both", direction="in", top=True, right=True)
    if label:
        ax.text(*loc, label, transform=ax.transAxes, va="top", ha="left", fontsize=8.5)


def headroom(ax, frac=0.18, log=False):
    """Extend the upper y limit so that the panel label does not touch the data."""
    lo, hi = ax.get_ylim()
    if log:
        ax.set_ylim(lo, hi * (hi / lo) ** frac)
    else:
        ax.set_ylim(lo, hi + frac * (hi - lo))


def save(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / f"{name}.pdf")
    plt.close(fig)
    print(f"wrote manuscript/figures/{name}.pdf")


def load():
    """Index the validation table by (quantity, parameters, beta, cos theta)."""
    table = {}
    for rec in records():
        key = (rec["quantity"], (rec["mA"], rec["muR"], rec["Nc"]), rec["beta"], rec["costheta"])
        table[key] = rec
    return table


def fig_validation(table):
    fig, ax = plt.subplots(2, 2, figsize=(WIDE, 4.6), sharex="col")
    plt.subplots_adjust(left=0.08, right=0.985, bottom=0.09, top=0.985, wspace=0.26, hspace=0.08)
    for col, quantity in enumerate(("M1_nf5", "M2_nf5")):
        worst = defaultdict(lambda: [0.0, 0.0])
        for (q, _, beta, _), rec in table.items():
            if q == quantity:
                for j, order in enumerate((4, 6)):
                    err = float(relative_error(rec[f"series{order}"], rec["exact"]))
                    worst[beta][j] = max(worst[beta][j], err)
        beta = np.array([float(b) for b in sorted(worst)])
        e4 = np.array([worst[b][0] for b in sorted(worst)])
        e6 = np.array([worst[b][1] for b in sorted(worst)])
        s4, s6 = (np.polyfit(np.log(beta), np.log(e), 1)[0] for e in (e4, e6))
        anchor = int(np.argmin(abs(beta - 0.02)))

        top = ax[0, col]
        top.loglog(beta, e4, color=BLACK, marker="o", mfc="none", label=r"through $\beta_A^4$")
        top.loglog(beta, e6, color=RED, marker="s", mfc="none", ls="--", label=r"through $\beta_A^6$")
        top.loglog(beta, e4[anchor] * (beta / beta[anchor]) ** 6, color=GREY, ls=":", lw=0.9, label=r"$\propto\beta_A^6$")
        top.loglog(beta, e6[anchor] * (beta / beta[anchor]) ** 8, color=GREY, ls="-.", lw=0.9, label=r"$\propto\beta_A^8$")
        top.set_ylabel("maximum relative error")
        headroom(top, 0.12, log=True)
        finish_axis(top, ("(a)" if col == 0 else "(b)") + rf"  $\mathcal{{M}}_{col + 1}$")
        top.text(0.04, 0.79, rf"fitted slopes ${s4:.3f}$, ${s6:.3f}$", transform=top.transAxes, fontsize=7.5)
        if col == 1:
            top.legend(loc="lower right")

        bottom = ax[1, col]
        bottom.semilogx(beta, e4 / beta**6, color=BLACK, marker="o", mfc="none", label=r"$E_4/\beta_A^6$")
        bottom.semilogx(beta, e6 / beta**8, color=RED, marker="s", mfc="none", ls="--", label=r"$E_6/\beta_A^8$")
        bottom.set_xlabel(r"$\beta_A$")
        bottom.set_ylabel("scaled relative error")
        headroom(bottom, 0.22)
        bottom.legend(loc="upper center")
        finish_axis(bottom, "(c)" if col == 0 else "(d)")
    save(fig, "fig1_validation")


def flavour_terms(table, i):
    """Exact n_f^0, 5 n_f^1 and 25 n_f^2 terms and the complete sum at the benchmark point."""
    betas = sorted({k[2] for k in table})
    rows = []
    for b in betas:
        key = (BENCHMARK, b, BENCHMARK_ANGLE)
        terms = [5**k * table[(f"M{i}_nf{k}",) + key]["exact"] for k in (0, 1, 2)]
        rows.append((float(b), terms, table[(f"M{i}_nf5",) + key]["exact"]))
    return rows


def fig_flavour(table):
    """Single-column figure: the two form factors stacked."""
    fig, axes = plt.subplots(2, 1, figsize=(3.375, 3.9), sharex=True)
    plt.subplots_adjust(left=0.19, right=0.97, bottom=0.11, top=0.985, hspace=0.07)
    styles = [(BLACK, "o", "-", r"$a_0=\mathcal{M}_{i,0}^{(2)}$"), (RED, "s", "--", r"$a_1=5\,\mathcal{M}_{i,1}^{(2)}$"),
              (BLUE, "^", "-.", r"$a_2=25\,\mathcal{M}_{i,2}^{(2)}$")]
    for j, i in enumerate((1, 2)):
        rows = flavour_terms(table, i)
        x = [r[0] for r in rows]
        for k, (c, m, ls, lab) in enumerate(styles):
            axes[j].loglog(x, [float(abs(r[1][k]) / abs(r[2])) for r in rows], color=c, marker=m, mfc="none", ls=ls, label=lab)
        axes[j].axhline(1.0, color=GREY, ls=":", lw=0.8)
        axes[j].set_ylabel(rf"$|a_r|/|\mathcal{{M}}_{i}(5)|$")
        headroom(axes[j], 0.3, log=True)
        finish_axis(axes[j], "(a)" if j == 0 else "(b)")
    axes[1].legend(loc="center left", bbox_to_anchor=(0.0, 0.40))
    axes[1].set_xlabel(r"$\beta_A$")
    save(fig, "fig2_flavour")


def fig_interference(table):
    fig, axes = plt.subplots(1, 2, figsize=(WIDE, 2.75), sharex=True, sharey=True)
    plt.subplots_adjust(left=0.085, right=0.985, bottom=0.18, top=0.76, wspace=0.16)
    for j, i in enumerate((1, 2)):
        x, deficit, interference = [], [], []
        for b, terms, total in flavour_terms(table, i):
            size = [abs(t) for t in terms]
            x.append(b)
            deficit.append(float(1 - abs(total) / sum(size)))
            interference.append(float((abs(total) ** 2 - sum(s**2 for s in size)) / sum(s**2 for s in size)))
        handles = [axes[j].semilogx(x, deficit, color=INK, marker="o", label=r"magnitude deficit $D_i$")[0],
                   axes[j].semilogx(x, interference, color=WINE, marker="s", ls="--",
                                    label=r"signed interference $I_i$")[0]]
        axes[j].axhline(0.0, color=GREY, lw=0.75)
        axes[j].set_ylim(-0.62, 0.55)
        axes[j].set_xlabel(r"$\beta_A$")
        if j == 0:
            axes[j].set_ylabel(r"$D_i$ or $I_i$")
        axes[j].set_title(rf"$\mathcal{{M}}_{i}$", pad=4)
        finish_axis(axes[j], "(a)" if j == 0 else "(b)")
    fig.legend(handles, [h.get_label() for h in handles], loc="upper center",
               bbox_to_anchor=(0.535, 0.98), ncol=2, frameon=False, handlelength=2.5)
    save(fig, "fig3_interference")


def fig_ratio(table):
    fig, axes = plt.subplots(1, 2, figsize=(WIDE, 2.85), sharex=True)
    plt.subplots_adjust(left=0.09, right=0.985, bottom=0.18, top=0.75, wspace=0.29)
    betas = sorted({k[2] for k in table})
    handles = []
    for angle, color, marker, ls in ((Fraction(0), INK, "o", "-"), (Fraction(1, 3), WINE, "s", "--"),
                                     (Fraction(2, 3), TEAL, "^", "-.")):
        ratios = [complex(table[("M2_nf5", BENCHMARK, b, angle)]["exact"] /
                          table[("M1_nf5", BENCHMARK, b, angle)]["exact"]) for b in betas]
        x = [float(b) for b in betas]
        phase = np.degrees(np.unwrap([math.atan2(r.imag, r.real) for r in ratios]))
        handles.append(axes[0].semilogx(x, [abs(r) for r in ratios], color=color, marker=marker, ls=ls,
                                        label=rf"$c={angle}$")[0])
        axes[1].semilogx(x, phase, color=color, marker=marker, ls=ls)
    axes[0].set_ylabel(r"$|\mathcal{M}_2/\mathcal{M}_1|$")
    axes[1].set_ylabel(r"$\arg(\mathcal{M}_2/\mathcal{M}_1)$ [deg]")
    for j, ax in enumerate(axes):
        ax.set_xlabel(r"$\beta_A$")
        finish_axis(ax, "(a)" if j == 0 else "(b)")
    fig.legend(handles, [h.get_label() for h in handles], loc="upper center",
               bbox_to_anchor=(0.535, 0.975), ncol=3, frameon=False, handlelength=2.4)
    save(fig, "fig4_ratio")


def fig_dps():
    """Fig. 4: (a) one- and two-loop forward coefficients against the helicity weights (Eqs. FA, FS);
    (b) their relative deviations from the prediction."""
    import json
    from matplotlib.lines import Line2D
    from ggaa import forward as F
    path = RESULTS / "dps_check.json"
    if not path.exists():
        print("skipped fig7_dps: run make dps first")
        return
    record = json.load(open(path))
    hh_path = RESULTS / "dps_hh.json"
    hh = json.load(open(hh_path))["points"] if hh_path.exists() else []
    rows = F.load() if F.POINTS.exists() else []
    colour = {(3, 5): RED, (4, 2): BLUE}
    fig, ax = plt.subplots(2, 1, figsize=(3.375, 4.3), sharex=True, gridspec_kw=dict(height_ratios=[1.35, 1]))
    plt.subplots_adjust(left=0.17, right=0.97, bottom=0.09, top=0.985, hspace=0.07)
    xs = np.linspace(0.0, 0.75, 200)
    x = (1 + xs) / 2
    ax[0].plot(xs, 1 / (x * (1 - x)) - 2 * x * (1 - x), color=BLACK, lw=1.0)
    ax[0].plot(xs, 1 / (x * (1 - x)) + 2 * x * (1 - x), color=BLACK, ls="--", lw=1.0)
    for (Nc, nf), c in colour.items():
        pts = [p for p in record["points"] if p["Nc"] == Nc and p["nf"] == nf]
        b = [p["beta"] for p in pts]
        d = [p["slope"] / (np.pi * Nc) for p in pts]
        ax[0].plot(b, d, "o", color=c, ms=5.5 if Nc == 3 else 3.5, mfc="none", mew=0.9, ls="none")
        ax[1].semilogy(b, [max(abs(di / p["F_A"] - 1), 1e-16) for di, p in zip(d, pts)], "o", color=c,
                       ms=5.5 if Nc == 3 else 3.5, mfc="none", mew=0.9, ls="none")
    if hh:
        ax[0].plot([p["beta"] for p in hh], [p["HH_tree_sign"] for p in hh], "^", color=GREEN, ms=5.0,
                   mfc="none", mew=0.9, ls="none")
        ax[1].semilogy([p["beta"] for p in hh], [max(abs(p["dev_HH"]), 1e-16) for p in hh], "^", color=GREEN,
                       ms=5.0, mfc="none", mew=0.9, ls="none")
    if rows:
        narrow, wide = F.table(rows, (7, 8, 9, 10)), F.table(rows, (6, 7, 8, 9, 10))
        for t, filled in ((narrow, True), (wide, False)):
            for (Nc, nf), c in colour.items():
                sel = [r for r in t if (r["Nc"], r["nf"]) == (Nc, nf)]
                b = [float(r["beta"]) + (-0.008 if Nc == 3 else 0.008) for r in sel]   # sets side by side
                for key, marker in (("3", "D"), ("2", "s")):
                    if filled:
                        fa = [float(F_A_value) for F_A_value in
                              (1 / (((1 + bb) / 2) * ((1 - bb) / 2)) - 2 * ((1 + bb) / 2) * ((1 - bb) / 2)
                               for bb in (float(r["beta"]) for r in sel))]
                        ax[0].plot(b, [r["ratio" + key] * f for r, f in zip(sel, fa)], marker, color=c,
                                   ms=3.2, mew=0.0, ls="none", alpha=0.85)
                    ax[1].semilogy(b, [r["dev" + key] for r in sel], marker, color=c, ms=3.4,
                                   mfc=c if filled else "none", mew=0.8, ls="none")
    ax[0].set_ylabel("normalised forward coefficient")
    ax[1].set_ylabel(r"$|\mathrm{fit}/\mathrm{prediction}-1|$")
    ax[1].set_xlabel(r"$\beta$")
    ax[1].set_ylim(3e-17, 3e-3)
    handles = [Line2D([], [], color=BLACK, lw=1.0, label=r"$F_A$"),
               Line2D([], [], color=BLACK, lw=1.0, ls="--", label=r"$F_S$"),
               Line2D([], [], color=GREY, marker="o", mfc="none", ls="none", label="one loop, $AA$"),
               Line2D([], [], color=GREEN, marker="^", mfc="none", ls="none", label="one loop, $HH$"),
               Line2D([], [], color=GREY, marker="D", ms=3.2, ls="none", label=r"two loops, $\mathcal{G}^3$"),
               Line2D([], [], color=GREY, marker="s", ms=3.2, ls="none", label=r"two loops, $\mathcal{G}^2$")]
    ax[0].legend(handles=handles, loc="upper left", fontsize=7, ncol=2, columnspacing=0.8, handlelength=1.6)
    finish_axis(ax[0], None)
    finish_axis(ax[1], None)
    ax[0].text(0.97, 0.06, "(a)", transform=ax[0].transAxes, ha="right", va="bottom", fontsize=8.5)
    ax[1].text(0.97, 0.06, "(b)", transform=ax[1].transAxes, ha="right", va="bottom", fontsize=8.5)
    save(fig, "fig7_dps")

def main():
    style()
    table = load()
    fig_validation(table)
    fig_flavour(table)
    fig_interference(table)
    fig_ratio(table)
    fig_dps()


if __name__ == "__main__":
    main()
