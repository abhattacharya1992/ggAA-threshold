#!/usr/bin/env python3
"""Why the earlier 2019 two-loop run is not an alternative to the published remainder (App. B).

Reads the 2019 files in data/provenance_2019 and checks:
1. the published renormalisation (Renorm1_N.m) uses Z_GG = Z_as exactly;
2. the earlier run's constant Z2 (archive two_loop/Renorm1.m) exceeds Z_GG^2 by exactly
   2 beta1 a_s^2/eps, both if its coupling variable is read as the bare coupling and if it is
   re-expanded in the renormalised one; the excess contains C_A^2, which no O_J term can supply;
3. the earlier run has no O_G-O_J counterterm; the published one has 64 Z_GG^2 and (8/6) Z_GG Z_GJ;
4. both runs use the same infrared subtraction (Catani2L1.m, Catani2L1_N.m): apart from the input
   file name, the test point, the form of zeta_3 and the output lines, the scripts coincide;
5. the June 2019 record (1byeps_Coef_Num): in both projections and at every kinematic point the
   C_F n_f single pole vanishes only for FLAG = 1/96, and 64 x 2 x (1/96) = 8/6, where 2 is the
   counterterm coefficient of Renorm2L_P1_old.m.
Writes results/renorm_2019.json; exits with an error if any check fails.
"""

import json
import pathlib
import re
import sys

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROV = ROOT / "data" / "provenance_2019"
TWO = PROV / "archive" / "two_loop"
RESULTS = ROOT / "results"

a, e, CA, CF, nf, ab = sp.symbols("a_s epsilon C_A C_F n_f a_b")
b0 = (11 * CA - 2 * nf) / 3
b1 = sp.Rational(34, 3) * CA**2 - 2 * nf * CF - sp.Rational(10, 3) * CA * nf


def mathematica(expr, **names):
    """Minimal Mathematica-to-sympy conversion for the polynomial expressions of these scripts."""
    s = re.sub(r"\(\*.*?\*\)", "", expr, flags=re.S)
    s = re.sub(r"(\d)\*\^(-?\d+)", r"\1e\2", s)                  # 1.2*^-75 -> 1.2e-75
    s = re.sub(r"\bas\b", "a_s", s.strip())                     # "as" is a Python keyword
    s = re.sub(r"(?<=[\w)])\s+(?=[A-Za-z(\d])", "*", s)          # implicit products "484 Ca^2", "Ca nf"
    s = s.replace("^", "**")
    table = {"Ca": CA, "Cf": CF, "nf": nf, "eps": e, "I": sp.I}
    table.update(names)
    return sp.sympify(s, locals=table)


def assignment(text, name):
    m = re.search(rf"^{name}\s*=\s*(.*?);", text, flags=re.S | re.M)
    if not m:
        raise SystemExit(f"{name} not found")
    return m.group(1)


def main():
    problems, out = [], {}
    new = (TWO / "Renorm1_N.m").read_text()
    old = (TWO / "Renorm1.m").read_text()

    # 1. published Z_GG = Z_as
    ZGG = mathematica(assignment(new, "ZGG"), a_s=a)
    Zas = mathematica(assignment(new, "Zas"), a_s=a, beta0=b0, beta1=b1)
    out["published_ZGG_equals_Zas"] = sp.simplify(sp.expand(ZGG - Zas)) == 0

    # 2. earlier Z2 against Z_GG^2, two readings of its coupling variable ascr
    x = sp.Symbol("x")
    Z2 = mathematica(assignment(old, "Z2"), ascr=x)
    Zas_x = 1 + a * 2 * b0 / e + a**2 * (4 * b0**2 / e**2 + b1 / e)
    d1 = sp.expand(sp.series(sp.expand(Z2.subs(x, a * Zas_x) - ZGG**2), a, 0, 3).removeO())
    a_of_ab = ab - 2 * b0 * ab**2 / e + (4 * b0**2 / e**2 - b1 / e) * ab**3
    zgg2_bare = sp.series(sp.expand((ZGG**2).subs(a, a_of_ab)), ab, 0, 3).removeO()
    d2 = sp.expand(sp.series(sp.expand(Z2.subs(x, ab) - zgg2_bare), ab, 0, 3).removeO())
    out["earlier_Z2_excess_is_2beta1_over_eps"] = (sp.simplify(d1 - 2 * b1 * a**2 / e) == 0
                                                    and sp.simplify(d2 - 2 * b1 * ab**2 / e) == 0)
    out["earlier_Z2_excess"] = str(sp.factor(d1))

    # 3. counterterms
    out["earlier_has_OG_OJ_counterterm"] = "ZGJ" in old
    m = re.search(r"\(\((\d+)\*\(ZGG\^2\)\*gg2AA\)\+\(\((\d+)/(\d+)\)\*ZGG\*ZGJ\*gg2AAj\)\)", new.replace(" ", ""))
    out["published_factors"] = [int(m.group(1)), f"{m.group(2)}/{m.group(3)}"] if m else None
    oldpub = (PROV / "Renorm2L_P1_old.m").read_text().replace(" ", "")
    m_old = re.search(r"\((\d+)\*ZGG\*ZGJ\*gg2AAj\)", oldpub)
    old_coef = int(m_old.group(1)) if m_old else None
    out["counterterm_coefficient_Renorm2L_P1_old"] = old_coef

    # 4. same infrared subtraction
    skip = re.compile(r"^(zeta\[3\]|zeta2\s*=|gg2AAR\s*=\s*Get|[xzYZ]\s*=|diff)")
    def body(path):
        lines = [re.sub(r"\s+", " ", l).strip() for l in path.read_text().splitlines()]
        return [l for l in lines if l and not skip.match(l)]
    out["same_infrared_subtraction"] = body(TWO / "Catani2L1.m") == body(TWO / "Catani2L1_N.m")

    # 5. the June 2019 record
    rec = (TWO / "1byeps_Coef_Num").read_text()
    FLAG, N, Nf, M2 = sp.symbols("FLAG N Nf M2")
    points = []
    for block in re.split(r"\*{20,}", rec):
        pt = re.search(r"x=([\d/]+);\s*z=([\d/]+);", block)
        if not pt:
            continue
        entry = {"x": pt.group(1), "z": pt.group(2)}
        for proj in ("A", "B"):
            m = re.search(rf"Proj {proj} -> 1/eps Coefficient = (.*?)\n\s*\n", block + "\n\n", flags=re.S)
            E = mathematica(m.group(1), FLAG=FLAG, N=N, Nf=Nf, M2=M2).subs({N: 3, M2: 1})
            nf_part = sp.expand(sp.diff(E, Nf))
            root = sp.solve(sp.Eq(nf_part, 0), FLAG)[0]
            size = abs(complex(nf_part.subs(FLAG, 0)))
            rest = abs(complex(E.subs(Nf, 0)))
            entry[f"proj{proj}_flag_root"] = [float(sp.re(root)), float(sp.im(root))]
            entry[f"proj{proj}_rel_dev_from_1_96"] = abs(complex(root) * 96 - 1)
            entry[f"proj{proj}_other_poles_relative"] = rest / size
        points.append(entry)
    out["pole_record_points"] = points
    out["pole_record_n_points"] = len(points)
    devs = [p[f"proj{q}_rel_dev_from_1_96"] for p in points for q in "AB"]
    out["pole_record_max_rel_dev"] = max(devs)
    out["published_coefficient_consistent"] = (out["published_factors"] == [64, "8/6"] and old_coef == 2
                                               and sp.Rational(64) * old_coef / 96 == sp.Rational(8, 6))

    if not out["published_ZGG_equals_Zas"]: problems.append("published Z_GG != Z_as")
    if not out["earlier_Z2_excess_is_2beta1_over_eps"]: problems.append("Z2 excess is not 2 beta1/eps")
    if out["earlier_has_OG_OJ_counterterm"]: problems.append("earlier run has an O_G-O_J counterterm")
    if not out["same_infrared_subtraction"]: problems.append("infrared subtractions differ")
    if len(points) != 7 or max(devs) > 1e-8: problems.append("FLAG root is not 1/96 at all points")
    if not out["published_coefficient_consistent"]: problems.append("64 x 2/96 != published coefficient")

    print(f"published Z_GG = Z_as: {out['published_ZGG_equals_Zas']}")
    print(f"earlier Z2 - Z_GG^2 = {out['earlier_Z2_excess']} (= 2 beta1 a_s^2/eps, both readings: "
          f"{out['earlier_Z2_excess_is_2beta1_over_eps']})")
    print(f"earlier run has an O_G-O_J counterterm: {out['earlier_has_OG_OJ_counterterm']}; "
          f"published factors {out['published_factors']}")
    print(f"same infrared subtraction: {out['same_infrared_subtraction']}")
    print(f"June 2019 record: {len(points)} points x 2 projections, FLAG root = 1/96 to "
          f"{out['pole_record_max_rel_dev']:.1e}; 64 x {old_coef}/96 = 8/6: {out['published_coefficient_consistent']}")
    RESULTS.mkdir(exist_ok=True)
    json.dump(out, open(RESULTS / "renorm_2019.json", "w"), indent=1)
    if problems:
        sys.exit("PROBLEMS: " + "; ".join(problems))


if __name__ == "__main__":
    main()
