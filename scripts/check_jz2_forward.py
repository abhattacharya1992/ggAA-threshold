#!/usr/bin/env python3
"""Forward limit of the two-loop |J_z| = 2 amplitude M_2 at fixed velocity.

For each beta the published two-loop P2 remainder (divided by 64) is evaluated at
cos(theta) = 1 - 10^-e, and M_2^(2)/m_A^2 = A u^2 + B u + C with u = ln(1-cos theta) is fitted.
a2 = 4A is the coefficient of l^2, l = ln(2 p_T/m_A).  It is compared with -4 N_c K1, where K1
is the forward one-loop M_2^(1)/m_A^2; at threshold a2 should approach -2 pi N_c^2.
The largest e of each beta is also evaluated with --check-extra more digits as a precision check.
Points are appended to results/jz2_forward_points.jsonl at once; a rerun skips them.
"""
import argparse, json, multiprocessing as mpc, os, pathlib, resource, sys, time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "python"))
import mpmath as mp
from ggaa.normalization import FILE_FACTOR
from ggaa.paths import INPUTS, RESULTS

POINTS = RESULTS / "jz2_forward_points.jsonl"
SUMMARY = RESULTS / "jz2_forward.txt"


def key(beta, e, Nc, nf, dps):
    return f"{beta}|{e}|{Nc}|{nf}|{dps}"


def evaluate(task):
    limit = int(task["memory_gb"] * 1024 ** 3)
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))   # a runaway worker fails alone
    os.nice(10)
    from ggaa.twoloop import remainder
    mp.mp.dps = task["dps"]
    beta, c = mp.mpf(task["beta"]), 1 - mp.mpf(10) ** (-task["e"])
    t0 = time.time()
    v = remainder(2, 2, beta, c, mA=task["mA"], muR=task["mA"], Nc=task["Nc"], nf=task["nf"])
    v = v / FILE_FACTOR / mp.mpf(task["mA"]) ** 2
    out = {k: task[k] for k in ("beta", "e", "Nc", "nf", "mA", "dps")}
    out.update(re=mp.nstr(v.real, task["dps"] - 10), im=mp.nstr(v.imag, task["dps"] - 10),
               seconds=round(time.time() - t0))
    return out


def load():
    done = {}
    if POINTS.exists():
        for line in POINTS.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                done[key(r["beta"], r["e"], r["Nc"], r["nf"], r["dps"])] = r
    return done


def value(r):
    return mp.mpc(mp.mpf(r["re"]), mp.mpf(r["im"]))


def one_loop_constant(beta, Nc, nf, mA):
    from ggaa.oneloop import OneLoop
    mp.mp.dps = 120
    v = OneLoop()(2, mp.mpf(beta), 1 - mp.mpf(10) ** -18, mA=mA, muR=mA, Nc=Nc, nf=nf)
    return v / mp.mpf(mA) ** 2


def fit(rows):
    mp.mp.dps = 60
    us = [-r["e"] * mp.log(10) for r in rows]
    vs = [value(r) for r in rows]
    M = mp.matrix([[u ** 2, u, 1] for u in us])
    coef = mp.lu_solve(M.T * M, M.T * mp.matrix(vs))
    res = max(abs(vs[i] - (coef[0] * us[i] ** 2 + coef[1] * us[i] + coef[2])) for i in range(len(rows)))
    return coef[0], coef[1], coef[2], res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--betas", nargs="+", default=["0.05", "0.1", "0.2", "0.3", "0.5"])
    ap.add_argument("--exponents", nargs="+", type=int, default=[4, 5, 6, 7])
    ap.add_argument("--dps", type=int, default=160)
    ap.add_argument("--check-extra", type=int, default=40)
    ap.add_argument("--Nc", type=int, default=3)
    ap.add_argument("--nf", type=int, default=5)
    ap.add_argument("--mA", type=int, default=125)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--memory-gb", type=float, default=3.0)
    a = ap.parse_args()
    if not (INPUTS / "Finite_2Loop_P2").exists():
        sys.exit("published two-loop files missing: run 'make inputs' first")
    ncpu = os.cpu_count() or 1
    if a.workers > max(1, ncpu - 2):
        sys.exit(f"--workers {a.workers} leaves fewer than two of {ncpu} cores free")

    tasks = []
    for b in a.betas:
        tasks += [dict(beta=b, e=e, dps=a.dps) for e in a.exponents]
        tasks.append(dict(beta=b, e=max(a.exponents), dps=a.dps + a.check_extra))
    for t in tasks:
        t.update(Nc=a.Nc, nf=a.nf, mA=a.mA, memory_gb=a.memory_gb)
    done = load()
    todo = [t for t in tasks if key(t["beta"], t["e"], a.Nc, a.nf, t["dps"]) not in done]
    print(f"{len(tasks)} points, {len(tasks) - len(todo)} done, {len(todo)} to go, {a.workers} workers", flush=True)
    RESULTS.mkdir(exist_ok=True)
    if todo:
        with mpc.get_context("spawn").Pool(a.workers, maxtasksperchild=1) as pool:
            for r in pool.imap_unordered(evaluate, todo):
                with POINTS.open("a") as f:
                    f.write(json.dumps(r) + "\n")
                done[key(r["beta"], r["e"], r["Nc"], r["nf"], r["dps"])] = r
                print(f"  beta={r['beta']} 1-c=1e-{r['e']} dps={r['dps']}: {r['re'][:18]} + {r['im'][:18]} i ({r['seconds']} s)", flush=True)

    lines = [f"N_c={a.Nc} n_f={a.nf} m_A=mu_R={a.mA}; threshold value of a2: -2 pi N_c^2 = {mp.nstr(-2 * mp.pi * a.Nc ** 2, 8)}",
             f"{'beta':>6} {'Re a2':>14} {'Im a2':>11} {'-4 N_c K1':>13} {'ratio':>11} {'fit resid.':>11} {'check digits':>13}"]
    for b in a.betas:
        rows = [done[key(b, e, a.Nc, a.nf, a.dps)] for e in a.exponents if key(b, e, a.Nc, a.nf, a.dps) in done]
        hi = done.get(key(b, max(a.exponents), a.Nc, a.nf, a.dps + a.check_extra))
        lo = done.get(key(b, max(a.exponents), a.Nc, a.nf, a.dps))
        mp.mp.dps = 60
        digits = "-"
        if hi and lo:
            d = abs(value(hi) - value(lo)) / abs(value(hi))
            digits = mp.nstr(-mp.log10(d), 3) if d > 0 else ">50"
        if len(rows) < 3:
            lines.append(f"{b:>6}  fewer than three points; check digits {digits}")
            continue
        A, B, C, res = fit(rows)
        pred = -4 * a.Nc * one_loop_constant(b, a.Nc, a.nf, a.mA).real
        mp.mp.dps = 60
        lines.append(f"{b:>6} {mp.nstr((4 * A).real, 9):>14} {mp.nstr((4 * A).imag, 3):>11} {mp.nstr(pred, 9):>13} "
                     f"{mp.nstr((4 * A).real / pred, 8):>11} {mp.nstr(res, 2):>11} {digits:>13}")
    SUMMARY.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
