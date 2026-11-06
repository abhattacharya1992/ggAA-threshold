#!/usr/bin/env python3
"""Compare the series rebuilt by reproduce_series.sh with data/construction of the package.

Each file data/construction/M{i}_nf{k}_beta_{n}.wl is compared with the coefficient of
beta^n of the rebuilt series, after both are brought to the same form by ggaa.wolfram.

    python3 compare_with_package.py --work WORK --package PACKAGE
"""

import argparse
import lzma
import pathlib
import re
import sys
import time

import sympy as sp
from sympy.parsing.mathematica import parse_mathematica


def text(path):
    raw = lzma.open(path, "rt").read() if path.suffix == ".xz" else path.read_text()
    return re.sub(r"\(\*.*?\*\)", "", raw, flags=re.S).replace("\\\n", "")


def split_top(s):
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
    pairs = (item.split("->", 1) for item in split_top(body[2:-2]))
    return {k.strip().strip('"'): v.strip() for k, v in pairs}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work", type=pathlib.Path, required=True)
    ap.add_argument("--package", type=pathlib.Path, required=True)
    args = ap.parse_args()
    sys.path.insert(0, str(args.package / "src" / "python"))
    from ggaa import wolfram as W

    results = args.work / "remaining" / "results"
    known = args.work / "remaining" / "known"
    construction = args.package / "data" / "construction"

    def parse(expr):
        expr = expr.replace(*W._PRINCIPAL_LOG)
        return parse_mathematica(re.sub(r"\bbeta\b", "bta", expr)).subs(sp.Symbol("Lb"), W.L)

    def jets(label):
        jet = association(text(results / f"{label}_series6.wl"))["Jet"].strip()
        return {int(k): parse(v) for k, v in (item.split("->", 1) for item in split_top(jet[2:-2]))}

    def p1_nf1():
        expr = parse(association(text(known / "P1_nf1_p6.wl"))["Expression"])
        bta = sp.Symbol("bta")
        poly = sp.Poly(sp.expand(expr.subs(sp.log(bta), W.L)), bta)
        return {m[0]: c for m, c in zip(poly.monoms(), poly.coeffs())}

    def nf2(label):
        body = text(known / f"{label}_p6.wl").strip()
        args_ = split_top(body[len("SeriesData["):-1])
        ell, ct = sp.symbols("ell ct")
        rules = {ct: W.costh, ell: sp.log(W.m2) - sp.log(W.mu2)}
        return {int(args_[3]) + j: parse(c).subs(rules) for j, c in enumerate(split_top(args_[2].strip()[1:-1]))}

    def power(path):
        tag = path.stem.rsplit("_beta_", 1)[1]
        return -int(tag[1:]) if tag.startswith("m") else int(tag)

    sources = {(1, 0): lambda: jets("P1_nf0"), (2, 0): lambda: jets("P2_nf0"), (2, 1): lambda: jets("P2_nf1"),
               (1, 1): p1_nf1, (1, 2): lambda: nf2("P1_nf2"), (2, 2): lambda: nf2("P2_nf2")}
    failures = []
    for (i, k), load in sources.items():
        start = time.time()
        rebuilt = load()
        files = {power(p): p for p in construction.glob(f"M{i}_nf{k}_beta_*.wl")}
        same_files = 0
        for n in sorted(set(files) | set(rebuilt)):
            ours = W.reduce(W.read(files[n])) if n in files else 0
            new = W.reduce(rebuilt[n]) if rebuilt.get(n, 0) != 0 else 0
            same = W.is_zero(sp.expand(ours - new))
            if not same:
                failures.append(f"M{i}_nf{k}_beta_{n}")
            same_files += same and n in files
        print(f"M{i} n_f^{k}: {same_files} of {len(files)} construction files reproduced "
              f"({time.time() - start:.0f} s)", flush=True)
    if failures:
        print("NOT REPRODUCED:", ", ".join(failures))
        sys.exit(1)
    print("All construction files reproduced.")


if __name__ == "__main__":
    main()
