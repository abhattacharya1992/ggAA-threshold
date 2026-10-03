#!/usr/bin/env python3
"""Forward limit of the two-loop equal-helicity amplitude M_1 at fixed velocity.

For each beta the released two-loop P1 remainder (divided by 64) is evaluated at
cos(theta) = 1 - 10^-e, and M_1^(2)/m_A^2 = b3 u^3 + b2 u^2 + b1 u + b0, u = ln(1 - cos theta),
is fitted.  In the Gram-determinant variable
    G = (1/2) ln[(m_A^4 - t u)/m_A^4] = u/2 + g0,   g0 = (1/2) ln[8 beta^2/(1 - beta^2)^2] + i pi/2,
the coefficients of G^3 and G^2 are A3 = 8 b3 and A2 = 4 b2 - 3 A3 g0.  They are compared with
the DPS prediction
    A3 = -(4/3) N_c C1,   A2 = -2 beta0 C1,   C1 = 2 pi N_c F_A(x),  x = (1 + beta)/2,
    F_A(x) = 1/[x(1-x)] - 2 x(1-x),
which the threshold series satisfy exactly through beta^6 (scripts/check_nlp_dps.py).
The largest e of each beta is also evaluated with --check-extra more digits as a precision check.
Points are appended to results/m1_forward_points.jsonl at once; a rerun skips them.
--source series replaces the exact remainder by the threshold series (a fast test of the
fit and of the conversion to G; only meaningful at small beta).
"""
import argparse, json, multiprocessing as mpc, os, pathlib, resource, sys, time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "python"))
import mpmath as mp
from ggaa.normalization import FILE_FACTOR
from ggaa.paths import INPUTS, RESULTS


def paths(source, Nc, nf):
    tag = "" if source == "exact" else "_series_test"
    return RESULTS / f"m1_forward_points{tag}.jsonl", RESULTS / f"m1_forward{tag}_Nc{Nc}_nf{nf}.txt"


def key(beta, e, Nc, nf, dps):
    return f"{beta}|{e}|{Nc}|{nf}|{dps}"


def evaluate(task):
    limit = int(task["memory_gb"] * 1024 ** 3)
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    os.nice(10)
    mp.mp.dps = task["dps"]
    beta, c = mp.mpf(task["beta"]), 1 - mp.mpf(10) ** (-task["e"])
    t0 = time.time()
    if task["source"] == "exact":
        from ggaa.twoloop import remainder
        v = remainder(1, 2, beta, c, mA=task["mA"], muR=task["mA"], Nc=task["Nc"], nf=task["nf"])
        v = v / FILE_FACTOR
    else:
        from ggaa.series import ThresholdSeries
        v = ThresholdSeries()(1, beta, c, mA=task["mA"], muR=task["mA"], Nc=task["Nc"], nf=task["nf"])
    v = v / mp.mpf(task["mA"]) ** 2
    out = {k: task[k] for k in ("beta", "e", "Nc", "nf", "mA", "dps", "source")}
    out.update(re=mp.nstr(v.real, task["dps"] - 10), im=mp.nstr(v.imag, task["dps"] - 10),
               seconds=round(time.time() - t0))
    return out


def load(points):
    done = {}
    if points.exists():
        for line in points.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                done[key(r["beta"], r["e"], r["Nc"], r["nf"], r["dps"])] = r
    return done


def value(r):
    return mp.mpc(mp.mpf(r["re"]), mp.mpf(r["im"]))


def fit(rows):
    """Least-squares cubic in u; returns (b3, b2, b1, b0, largest residual)."""
    mp.mp.dps = 80
    us = [-r["e"] * mp.log(10) for r in rows]
    vs = [value(r) for r in rows]
    M = mp.matrix([[u ** 3, u ** 2, u, 1] for u in us])
    coef = mp.lu_solve(M.T * M, M.T * mp.matrix(vs))
    res = max(abs(vs[i] - sum(coef[k] * us[i] ** (3 - k) for k in range(4))) for i in range(len(rows)))
    return coef[0], coef[1], coef[2], coef[3], res


def prediction(beta, Nc, nf):
    mp.mp.dps = 80
    b = mp.mpf(beta)
    x = (1 + b) / 2
    C1 = 2 * mp.pi * Nc * (1 / (x * (1 - x)) - 2 * x * (1 - x))
    beta0 = mp.mpf(11 * Nc - 2 * nf) / 3
    return -mp.mpf(4) / 3 * Nc * C1, -2 * beta0 * C1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--betas", nargs="+", default=["0.05", "0.2", "0.3", "0.5", "0.7"])
    ap.add_argument("--exponents", nargs="+", type=int, default=[7, 8, 9, 10, 11])
    ap.add_argument("--dps", type=int, default=220)
    ap.add_argument("--check-extra", type=int, default=40)
    ap.add_argument("--Nc", type=int, default=3)
    ap.add_argument("--nf", type=int, default=5)
    ap.add_argument("--mA", type=int, default=125)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--memory-gb", type=float, default=3.0)
    ap.add_argument("--source", choices=["exact", "series"], default="exact")
    a = ap.parse_args()
    if a.source == "exact" and not (INPUTS / "Finite_2Loop_P1").exists():
        sys.exit("published two-loop files missing: run 'make inputs' first")
    ncpu = os.cpu_count() or 1
    if a.workers > max(1, ncpu - 2):
        sys.exit(f"--workers {a.workers} leaves fewer than two of {ncpu} cores free")
    points, summary = paths(a.source, a.Nc, a.nf)

    tasks = []
    for b in a.betas:
        tasks += [dict(beta=b, e=e, dps=a.dps) for e in a.exponents]
        tasks.append(dict(beta=b, e=max(a.exponents), dps=a.dps + a.check_extra))
    for t in tasks:
        t.update(Nc=a.Nc, nf=a.nf, mA=a.mA, memory_gb=a.memory_gb, source=a.source)
    done = load(points)
    todo = [t for t in tasks if key(t["beta"], t["e"], a.Nc, a.nf, t["dps"]) not in done]
    print(f"{len(tasks)} points, {len(tasks) - len(todo)} done, {len(todo)} to go, {a.workers} workers", flush=True)
    RESULTS.mkdir(exist_ok=True)
    if todo:
        with mpc.get_context("spawn").Pool(a.workers, maxtasksperchild=1) as pool:
            for r in pool.imap_unordered(evaluate, todo):
                with points.open("a") as f:
                    f.write(json.dumps(r) + "\n")
                done[key(r["beta"], r["e"], r["Nc"], r["nf"], r["dps"])] = r
                print(f"  beta={r['beta']} 1-c=1e-{r['e']} dps={r['dps']}: "
                      f"{r['re'][:18]} + {r['im'][:18]} i ({r['seconds']} s)", flush=True)

    lines = [f"source={a.source}  N_c={a.Nc} n_f={a.nf} m_A=mu_R={a.mA}",
             "A3, A2: coefficients of G^3, G^2 in M_1^(2)/m_A^2; prediction -(4/3)N_c C1 and -2 beta0 C1",
             f"{'beta':>6} {'Re A3':>14} {'pred A3':>14} {'Re/pred-1':>10} {'Im/|A3|':>9} "
             f"{'Re A2':>14} {'pred A2':>14} {'Re/pred-1':>10} {'Im/|A2|':>9} {'fit res.':>9} {'chk dig.':>8}"]
    for b in a.betas:
        rows = [done[key(b, e, a.Nc, a.nf, a.dps)] for e in a.exponents if key(b, e, a.Nc, a.nf, a.dps) in done]
        hi = done.get(key(b, max(a.exponents), a.Nc, a.nf, a.dps + a.check_extra))
        lo = done.get(key(b, max(a.exponents), a.Nc, a.nf, a.dps))
        mp.mp.dps = 80
        digits = "-"
        if hi and lo:
            d = abs(value(hi) - value(lo)) / abs(value(hi))
            digits = mp.nstr(-mp.log10(d), 3) if d > 0 else ">50"
        if len(rows) < 4:
            lines.append(f"{b:>6}  fewer than four points; check digits {digits}")
            continue
        b3, b2, b1, b0, res = fit(rows)
        bb = mp.mpf(b)
        g0 = mp.log(8 * bb ** 2 / (1 - bb ** 2) ** 2) / 2 + 1j * mp.pi / 2
        A3 = 8 * b3
        A2 = 4 * b2 - 3 * A3 * g0
        p3, p2 = prediction(b, a.Nc, a.nf)
        lines.append(f"{b:>6} {mp.nstr(A3.real, 10):>14} {mp.nstr(p3, 10):>14} {mp.nstr(A3.real / p3 - 1, 2):>10} "
                     f"{mp.nstr(A3.imag / abs(A3), 2):>9} {mp.nstr(A2.real, 10):>14} {mp.nstr(p2, 10):>14} "
                     f"{mp.nstr(A2.real / p2 - 1, 2):>10} {mp.nstr(A2.imag / abs(A2), 2):>9} "
                     f"{mp.nstr(res, 2):>9} {digits:>8}")
    summary.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
