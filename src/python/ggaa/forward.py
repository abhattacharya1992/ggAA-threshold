"""Forward limit c -> 1 at fixed velocity of the two-loop equal-helicity amplitude M_1.

The points are written by scripts/check_m1_forward.py to results/m1_forward_points.jsonl:
M_1^(2)/m_A^2 (package normalisation) at cos(theta) = 1 - 10^-e.  A cubic polynomial in
u = ln(1 - c) is fitted and rewritten in the Gram-determinant variable
    G = (1/2) ln[(m_A^4 - t u)/m_A^4] = u/2 + g0,  g0 = (1/2) ln[8 beta^2/(1 - beta^2)^2] + i pi/2,
so that the coefficients of G^3 and G^2 are A3 = 8 b3 and A2 = 4 b2 - 3 A3 g0.  They are compared
with C1 [-(4/3) N_c] and C1 [-2 beta_0], where C1 = 2 pi N_c F_A(x), x = (1 + beta)/2, is the
coefficient of G in the one-loop remainder.

The forward limit of the one-loop |J_z| = 2 remainder has the closed form
    K1 = pi N_c [x^4 + (1 - x)^4] / [x (1 - x)]   (units of m_A^2).
"""

import json

import mpmath as mp

from .paths import RESULTS

POINTS = RESULTS / "m1_forward_points.jsonl"
MAIN_DPS = 180
CHECK_DPS = 220


def load(path=POINTS):
    rows = []
    for line in path.read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            if r.get("source", "exact") == "exact":
                rows.append(r)
    return rows


def sets(rows):
    return sorted({(r["Nc"], r["nf"]) for r in rows})


def betas(rows, Nc, nf):
    return sorted({r["beta"] for r in rows if (r["Nc"], r["nf"]) == (Nc, nf)}, key=float)


def value(r):
    return mp.mpc(mp.mpf(r["re"]), mp.mpf(r["im"]))


def F_A(x):
    return 1 / (x * (1 - x)) - 2 * x * (1 - x)


def prediction(beta, Nc, nf):
    """Predicted coefficients of G^3 and G^2."""
    x = (1 + mp.mpf(beta)) / 2
    C1 = 2 * mp.pi * Nc * F_A(x)
    beta0 = mp.mpf(11 * Nc - 2 * nf) / 3
    return -mp.mpf(4) / 3 * Nc * C1, -2 * beta0 * C1


def k1_closed(beta, Nc):
    x = (1 + mp.mpf(beta)) / 2
    return mp.pi * Nc * (x**4 + (1 - x)**4) / (x * (1 - x))


def fit(rows, beta, Nc, nf, exponents, dps=MAIN_DPS):
    """Return (A3, A2) for one velocity and parameter set, fitted over the given exponents."""
    with mp.workdps(80):
        pick = {r["e"]: r for r in rows
                if (r["beta"], r["Nc"], r["nf"], r["dps"]) == (beta, Nc, nf, dps) and r["e"] in exponents}
        es = sorted(pick)
        if len(es) < 4:
            raise ValueError(f"fewer than four points for beta={beta}, Nc={Nc}, nf={nf}")
        us = [-e * mp.log(10) for e in es]
        vs = [value(pick[e]) for e in es]
        M = mp.matrix([[u**3, u**2, u, 1] for u in us])
        b3, b2 = mp.lu_solve(M.T * M, M.T * mp.matrix(vs))[:2]
        bb = mp.mpf(beta)
        g0 = mp.log(8 * bb**2 / (1 - bb**2) ** 2) / 2 + 1j * mp.pi / 2
        A3 = 8 * b3
        return A3, 4 * b2 - 3 * A3 * g0


def table(rows, exponents):
    """Relative deviations from the prediction for every set and velocity."""
    out = []
    for Nc, nf in sets(rows):
        for b in betas(rows, Nc, nf):
            A3, A2 = fit(rows, b, Nc, nf, exponents)
            p3, p2 = prediction(b, Nc, nf)
            out.append(dict(Nc=Nc, nf=nf, beta=b,
                            ratio3=float(A3.real / p3), ratio2=float(A2.real / p2),
                            dev3=float(abs(A3.real / p3 - 1)), dev2=float(abs(A2.real / p2 - 1)),
                            imag3=float(abs(A3.imag) / abs(A3)), imag2=float(abs(A2.imag) / abs(A2))))
    return out


def check_digits(rows):
    """Smallest number of agreeing digits between the two working precisions."""
    worst = 1e9
    with mp.workdps(CHECK_DPS + 10):
        for hi in (r for r in rows if r["dps"] == CHECK_DPS):
            for lo in (r for r in rows if r["dps"] == MAIN_DPS and
                       (r["beta"], r["e"], r["Nc"], r["nf"]) == (hi["beta"], hi["e"], hi["Nc"], hi["nf"])):
                d = abs(value(hi) - value(lo)) / abs(value(hi))
                worst = min(worst, float(-mp.log10(d)) if d > 0 else float(CHECK_DPS))
    return worst
