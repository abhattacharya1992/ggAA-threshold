#!/usr/bin/env python3
"""Check that data/construction and data/validation are the outputs of the construction runs.

data/runs holds the outputs of the construction and validation runs (xz-compressed), made
with the code in src/construction.  This script checks

1. the outputs against data/runs/SHA256SUMS (checksums of the uncompressed files);
2. every file data/construction/M{i}_nf{k}_beta_{n}.wl against the coefficient of beta^n in
   the corresponding run output: P1 n_f^0, P2 n_f^0, P2 n_f^1 from the recovery run
   (P*_series6.wl), P1 n_f^1 and the n_f^2 coefficients from the earlier component runs
   (known/*_p6.wl).  Both sides are brought to the same form by ggaa.wolfram.reduce
   and their difference must vanish identically;  powers that are absent on one side must
   vanish on the other;
3. data/validation/validation_1080.csv against the validation run (validation_all.csv):
   same 1080 records, and the six stored numbers agree to the 40 digits kept in the table.

Writes results/provenance.json.  Takes about two minutes (make provenance).
"""

import csv
import hashlib
import io
import json
import lzma
import pathlib
import re
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "python"))

import mpmath as mp
import sympy as sp
from sympy.parsing.mathematica import parse_mathematica

from ggaa import wolfram as W
from ggaa.paths import CONSTRUCTION, RESULTS

OUT = ROOT / "data" / "runs"
VALIDATION = ROOT / "data" / "validation" / "validation_1080.csv"


def check_integrity():
    bad = []
    lines = [l for l in (OUT / "SHA256SUMS").read_text().split("\n") if l]
    for line in lines:
        digest, name = line.split(None, 1)
        if hashlib.sha256(lzma.open(OUT / name).read()).hexdigest() != digest:
            bad.append(name)
    return len(lines), bad


def text(path):
    raw = lzma.open(path, "rt").read() if path.suffix == ".xz" else path.read_text()
    return re.sub(r"\(\*.*?\*\)", "", raw, flags=re.S).replace("\\\n", "")


def split_top(s):
    """Split at commas that are not nested in brackets, associations or strings."""
    out, cur, depth, i, in_string = [], [], 0, 0, False
    while i < len(s):
        c = s[i]
        if in_string:
            cur.append(c)
            in_string = not (c == '"' and s[i - 1] != "\\")
        elif c == '"':
            in_string = True
            cur.append(c)
        elif s.startswith("<|", i) or s.startswith("|>", i):
            depth += 1 if s[i] == "<" else -1
            cur.append(s[i:i + 2])
            i += 2
            continue
        elif c in "([{":
            depth += 1
            cur.append(c)
        elif c in ")]}":
            depth -= 1
            cur.append(c)
        elif c == "," and depth == 0:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(c)
        i += 1
    out.append("".join(cur))
    return out


def association(body):
    body = body.strip()
    assert body.startswith("<|") and body.endswith("|>")
    pairs = (item.split("->", 1) for item in split_top(body[2:-2]))
    return {k.strip().strip('"'): v.strip() for k, v in pairs}


def parse(expr):
    expr = expr.replace(*W._PRINCIPAL_LOG)
    expr = re.sub(r"\bbeta\b", "bta", expr)
    return parse_mathematica(expr).subs(sp.Symbol("Lb"), W.L)


def recovery_coefficients(label):
    jet = association(text(OUT / f"{label}_series6.wl.xz"))["Jet"].strip()
    return {int(k): parse(v) for k, v in (item.split("->", 1) for item in split_top(jet[2:-2]))}


def p1_nf1_coefficients():
    expr = parse(association(text(OUT / "known" / "P1_nf1_p6.wl.xz"))["Expression"])
    bta = sp.Symbol("bta")
    poly = sp.Poly(sp.expand(expr.subs(sp.log(bta), W.L)), bta)
    return {m[0]: c for m, c in zip(poly.monoms(), poly.coeffs())}


def nf2_coefficients(label):
    body = text(OUT / "known" / f"{label}_p6.wl.xz").strip()
    args = split_top(body[len("SeriesData["):-1])
    assert int(args[5]) == 1
    ell, ct = sp.symbols("ell ct")
    rules = {ct: W.costh, ell: sp.log(W.m2) - sp.log(W.mu2)}
    return {int(args[3]) + j: parse(c).subs(rules) for j, c in enumerate(split_top(args[2].strip()[1:-1]))}


def power(path):
    tag = path.stem.rsplit("_beta_", 1)[1]
    return -int(tag[1:]) if tag.startswith("m") else int(tag)


def compare_construction():
    sources = {(1, 0): lambda: recovery_coefficients("P1_nf0"),
               (2, 0): lambda: recovery_coefficients("P2_nf0"),
               (2, 1): lambda: recovery_coefficients("P2_nf1"),
               (1, 1): p1_nf1_coefficients,
               (1, 2): lambda: nf2_coefficients("P1_nf2"),
               (2, 2): lambda: nf2_coefficients("P2_nf2")}
    record, failures = {}, []
    for (i, k), load in sources.items():
        start = time.time()
        reference = load()
        files = {power(p): p for p in CONSTRUCTION.glob(f"M{i}_nf{k}_beta_*.wl")}
        checked = 0
        for n in sorted(set(files) | set(reference)):
            ours = W.reduce(W.read(files[n])) if n in files else 0
            theirs = W.reduce(reference[n]) if reference.get(n, 0) != 0 else 0
            same = W.is_zero(sp.expand(ours - theirs))
            if not same:
                failures.append(f"M{i}_nf{k}_beta_{n}")
            checked += same and n in files
        record[f"M{i}_nf{k}"] = {"files": len(files), "files_checked": checked}
        print(f"   M{i} n_f^{k}: {checked} of {len(files)} construction files identical to the run output"
              f" ({time.time() - start:.0f} s)", flush=True)
    return record, failures


def compare_validation():
    mp.mp.dps = 60
    run = list(csv.DictReader(io.StringIO(lzma.open(OUT / "validation_all.csv.xz", "rt").read())))
    ours = list(csv.DictReader(open(VALIDATION)))

    def run_key(r):
        m2, mu2, nc = (int(v) for v in r["Parameters"].strip("{}").split(","))
        name = r["Quantity"].replace("P1", "M1").replace("P2", "M2").replace("_FULL_nf5", "_nf5")
        return name, int(mp.sqrt(m2)), int(mp.sqrt(mu2)), nc, r["Beta"], r["CosTheta"]

    def our_key(r):
        return (r["quantity"], int(r["mA_GeV"]), int(r["muR_GeV"]), int(r["Nc"]), r["beta"], r["costheta"])

    run = {run_key(r): r for r in run}
    ours = {our_key(r): r for r in ours}
    columns = [("re_exact", "ReUnexpanded"), ("im_exact", "ImUnexpanded"), ("re_series4", "ReSeries4"),
               ("im_series4", "ImSeries4"), ("re_series6", "ReSeries6"), ("im_series6", "ImSeries6")]
    worst = mp.mpf(0)
    for key in set(run) & set(ours):
        for a, b in columns:
            x, y = mp.mpf(ours[key][a]), mp.mpf(run[key][b])
            if x != 0:
                worst = max(worst, abs(x - y) / abs(x))
    return {"records": len(ours), "matched": len(set(run) & set(ours)),
            "only_in_table": len(set(ours) - set(run)), "only_in_run": len(set(run) - set(ours)),
            "max_rel_difference": mp.nstr(worst, 3)}


def main():
    print("1. integrity of data/runs")
    count, bad = check_integrity()
    print(f"   {count} files, {len(bad)} with a wrong checksum {bad if bad else ''}")
    print("2. data/construction against the run outputs")
    construction, failures = compare_construction()
    print("3. data/validation against the validation run")
    validation = compare_validation()
    print(f"   {validation['matched']} of {validation['records']} records matched; largest relative "
          f"difference {validation['max_rel_difference']} (the table keeps 40 digits)")
    ok = (not bad and not failures and validation["matched"] == validation["records"] == 1080
          and validation["only_in_run"] == 0 and mp.mpf(validation["max_rel_difference"]) < mp.mpf("1e-38"))
    RESULTS.mkdir(exist_ok=True)
    json.dump({"integrity": {"files": count, "bad": bad}, "construction": construction,
               "construction_failures": failures, "validation": validation, "passed": ok},
              open(RESULTS / "provenance.json", "w"), indent=1)
    print("All checks passed." if ok else f"FAILED: {failures or bad or validation}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
