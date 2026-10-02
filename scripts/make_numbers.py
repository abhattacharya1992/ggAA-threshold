#!/usr/bin/env python3
"""Write manuscript/numbers.tex: every number quoted in the text, computed from the data
folders and the results files.  Nothing in numbers.tex is entered by hand.

Sources:
  data/validation/validation_1080.csv, forward_angles.csv   grid, Bose symmetry, phases, flavour sums
  results/accuracy.json          maximum errors and slopes       (make validation)
  results/leading_logs.json      leading-power coefficients      (make leading-logs)
  results/one_loop.json          one-loop threshold behaviour    (make one-loop)
  results/mb_constants.json      PSLQ reductions                 (make constants)
  results/normalization.json     normalization checks            (make normalization)
  results/exact_amplitude.txt    independent evaluation          (make exact-check)
  results/lo_cross_section.csv, virtual_correction.csv           (make lo, make virtual)
  data/timing/                   timing run
  data/series, inputs/           file sizes
  src/python/ggaa/benchmark.py, normalization.py                 parameters
"""

import csv
import json
import math
import pathlib
import re
import sys
from collections import defaultdict
from fractions import Fraction

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import mpmath as mp
import sympy as sp

from ggaa import benchmark
from ggaa.normalization import FILE_FACTOR
from ggaa.paths import DATA, INPUTS, RESULTS, SERIES
from ggaa.validation import records

ROOT = DATA.parent
WORDS = "zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen".split()
macros = {}


def define(name, value):
    macros[name] = str(value)


def sci(x, digits=3):
    """\\ensuremath{m\\times10^{e}} with `digits` significant digits."""
    if x == 0:
        return "0"
    exponent = math.floor(math.log10(abs(x)))
    mantissa = x / 10**exponent
    text = f"{mantissa:.{digits - 1}f}"
    if text.startswith("10"):
        exponent += 1
        text = f"{mantissa / 10:.{digits - 1}f}"
    return f"\\ensuremath{{{text}\\times10^{{{exponent}}}}}"


def fixed(x, digits):
    return f"{x:.{digits}f}"


def frac(x):
    f = Fraction(x).limit_denominator(10**6)
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


def power_of_ten(x):
    e = round(math.log10(x))
    return f"10^{{{e}}}" if abs(10**e - x) < 1e-12 * x else frac(x)


def listing(items):
    items = list(items)
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def rational_times(value, unit, tolerance=1e-8):
    """Recognize value = (p/q) * unit; return (Fraction, relative deviation)."""
    f = Fraction(value / unit).limit_denominator(100)
    return f, abs(value - float(f) * unit) / abs(value)


def tex_fraction(f, symbol):
    sign = "-" if f < 0 else ""
    f = abs(f)
    if f.denominator == 1:
        return f"{sign}{f.numerator if f.numerator != 1 else ''}{symbol}"
    return f"{sign}\\tfrac{{{f.numerator}}}{{{f.denominator}}}{symbol}"


# ---------------------------------------------------------------- benchmark parameters
RECOMMENDED_BETA = Fraction(1, 5)       # upper end of the recommended range
EDGE_BETA = Fraction(3, 10)             # largest validated velocity


def parameters():
    define("RecommendedBeta", f"{float(RECOMMENDED_BETA):g}")
    define("EdgeBeta", f"{float(EDGE_BETA):g}")
    define("BenchMA", benchmark.M_A)
    define("BenchMuR", benchmark.MU_R)
    define("BenchNc", benchmark.N_C)
    define("BenchNf", benchmark.N_F)
    define("BenchCos", frac(benchmark.COS_THETA))
    define("BenchAlphaS", benchmark.ALPHA_S)
    define("BenchAs", f"{benchmark.ALPHA_S / (4 * math.pi):.4f}")
    define("BenchCotBeta", benchmark.COT_BETA)
    define("FileFactor", FILE_FACTOR)
    define("FileFactorRoot", int(round(math.sqrt(FILE_FACTOR))))


# ---------------------------------------------------------------- validation grid and accuracy
def validation():
    rows = list(records())
    betas = sorted({r["beta"] for r in rows})
    angles = sorted({r["costheta"] for r in rows})
    params = {(r["mA"], r["muR"], r["Nc"]) for r in rows}
    define("ValidationRows", len(rows))
    define("BetaPointCount", len(betas))
    define("AnglePointCount", len(angles))
    define("ParameterSetCount", len(params))
    define("GridPointCount", len(betas) * len(angles) * len(params))
    define("ValBetaMin", f"\\ensuremath{{{power_of_ten(float(betas[0]))}}}")
    define("ValBetaMax", frac(betas[-1]))
    define("ValCosMax", frac(max(angles)))

    acc = json.load(open(RESULTS / "accuracy.json"))
    err = acc["max_error"]
    for i, name in ((1, "Pone"), (2, "Ptwo")):
        e6 = {Fraction(b): v for b, v in err[f"M{i}_nf5_order6"].items()}
        define(f"{name}ErrLePointTwo", sci(max(v for b, v in e6.items() if b <= RECOMMENDED_BETA)))
        define(f"{name}ErrAtPointThree", sci(e6[EDGE_BETA]))
        define(f"{name}SlopeFour", fixed(acc["slopes"][f"M{i}_nf5_order4"], 3))
        define(f"{name}SlopeSix", fixed(acc["slopes"][f"M{i}_nf5_order6"], 3))

    # Bose symmetry: c against -c for the exact values and both truncations
    table = {(r["quantity"], r["mA"], r["muR"], r["Nc"], r["beta"], r["costheta"]): r for r in rows}
    exact_checks, series_checks, worst = 0, 0, 0
    for key, r in table.items():
        if r["costheta"] > 0:
            partner = table[key[:5] + (-r["costheta"],)]
            exact_checks += 1
            worst = max(worst, float(abs(r["exact"] - partner["exact"]) / abs(r["exact"])))
            for o in (4, 6):
                series_checks += 1
                worst = max(worst, float(abs(r[f"series{o}"] - partner[f"series{o}"]) / abs(r[f"series{o}"])))
    define("AngularReferenceChecks", exact_checks)
    define("AngularSeriesChecks", series_checks)
    define("BoseStatement", "are identical in all stored digits" if worst == 0 else f"agree to {sci(worst, 2)}")

    # relative magnitude and phase of M2/M1 at the benchmark parameters
    bench = (benchmark.M_A, benchmark.MU_R, benchmark.N_C)
    c0 = benchmark.COS_THETA
    ratio = lambda b, c: table[("M2_nf5",) + bench + (b, c)]["exact"] / table[("M1_nf5",) + bench + (b, c)]["exact"]
    lo, hi = Fraction(1, 10), Fraction(3, 10)
    define("RatioBetaLow", frac(lo))
    define("RatioBetaHigh", frac(hi))
    define("RhoLow", fixed(float(abs(ratio(lo, c0))), 4))
    define("RhoHigh", fixed(float(abs(ratio(hi, c0))), 4))
    define("PhaseLow", fixed(float(mp.degrees(mp.arg(ratio(lo, c0)))), 2))
    define("PhaseHigh", fixed(float(mp.degrees(mp.arg(ratio(hi, c0)))), 2))
    spreads = {b: max(float(mp.degrees(mp.arg(ratio(b, c)))) for c in angles) -
                  min(float(mp.degrees(mp.arg(ratio(b, c)))) for c in angles) for b in betas}
    b_max = max(spreads, key=spreads.get)
    define("PhaseSpread", fixed(spreads[b_max], 2))
    define("PhaseSpreadBeta", frac(b_max))
    for i, name in ((1, "One"), (2, "Two")):
        for b, tag in ((RECOMMENDED_BETA, "PointTwo"), (EDGE_BETA, "PointThree")):
            worst_phase = max(float(abs(mp.degrees(mp.arg(r["series6"] / r["exact"]))))
                              for r in rows if r["quantity"] == f"M{i}_nf5" and r["beta"] == b)
            define(f"PhaseErr{name}{tag}", sci(worst_phase, 3))

    # flavour decomposition at the benchmark point
    deficits, inter = defaultdict(list), defaultdict(list)
    for b in betas:
        for i in (1, 2):
            terms = [5**k * table[(f"M{i}_nf{k}",) + bench + (b, c0)]["exact"] for k in (0, 1, 2)]
            total = sum(terms)
            size = [abs(t) for t in terms]
            deficits[i].append(float(1 - abs(total) / sum(size)))
            inter[i].append((b, float((abs(total) ** 2 - sum(s**2 for s in size)) / sum(s**2 for s in size))))
    define("PoneCancellationMin", f"{100 * min(deficits[1]):.1f}\\%")
    define("PoneCancellationMax", f"{100 * max(deficits[1]):.1f}\\%")
    define("ItwoPositiveBetas", listing([f"${frac(b)}$" for b, v in inter[2] if v > 0]))
    define("ItwoNegativeFrom", frac(min(b for b, v in inter[2] if v < 0)))
    assert all(v < 0 for b, v in inter[1]), "I_1 changes sign: revise the text of Sec. 5.3"
    assert max(b for b, v in inter[2] if v > 0) < min(b for b, v in inter[2] if v < 0), "I_2 sign pattern changed"
    return err


# ---------------------------------------------------------------- forward angles
def forward(err):
    rows = list(csv.DictReader(open(DATA / "validation" / "forward_angles.csv")))
    betas = sorted({Fraction(r["beta"]) for r in rows})
    cosines = sorted({Fraction(r["costheta"]) for r in rows}, key=lambda c: (c < 0, c))
    define("ForwardBetaList", listing([f"${float(b):g}$" for b in betas]))
    define("ForwardCosList", listing([f"${float(c):g}$" for c in cosines]))
    define("ForwardCosMax", f"{float(max(cosines)):g}")
    exceptions = []
    for q in ("M1", "M2"):
        for b in betas:
            grid = err[f"{q}_nf5_order6"][str(b)]
            worst = max((float(r["relative_error6"]), r["costheta"]) for r in rows
                        if r["quantity"] == f"{q}_nf5" and Fraction(r["beta"]) == b)
            if worst[0] > grid:
                exceptions.append(f"$\\mathcal M_{q[1]}$ at $\\bA={float(b):g}$, where the largest value, "
                                  f"{sci(worst[0])} at $c={float(Fraction(worst[1])):g}$, exceeds it by "
                                  f"{100 * (worst[0] / grid - 1):.0f}\\%")
    define("ForwardSummary", "All values are below the maximum on the main grid at the same $\\bA$"
           + (", except for " + listing(exceptions) if exceptions else "") + ".")


# ---------------------------------------------------------------- independent evaluation
def exact_check():
    text = (RESULTS / "exact_amplitude.txt").read_text()
    full = re.findall(r"(M\d)\s+\(mA, muR, Nc\) = \((\d+), (\d+), (\d+)\)\s+beta =\s+(\S+)\s+cos =\s+(\S+):\s+stable at (\d+) digits, "
                      r"relative difference\s+(\S+)", text)
    comps = re.findall(r"M\d, n_f\^\d: relative difference (\S+)", text)
    head = re.search(r"Flavour coefficients at \(mA, muR, Nc\) = \((\d+), (\d+), (\d+)\), beta = (\S+), cos = (\S+) \((\d+) digits\)", text)
    points = {row[1:6] for row in full}
    define("ExactPointCount", WORDS[len(points)])
    define("ExactMaxComplete", sci(max(float(r[7]) for r in full), 2))
    define("ExactMaxComponent", sci(max(float(c) for c in comps), 2))
    define("ExactComponentCount", WORDS[len(comps)])
    define("ExactComponentBeta", head.group(4).rstrip(","))
    define("ExactComponentCos", head.group(5))
    define("ExactDigits", max(int(r[6]) for r in full))
    define("ExactParameterSets", WORDS[len({row[1:4] for row in full})])


# ---------------------------------------------------------------- one loop
def one_loop():
    d = json.load(open(RESULTS / "one_loop.json"))
    define("OneLoopFitBetas", listing([f"$10^{{{int(round(math.log10(float(b))))}}}$" for b in d["fit_betas"]]))
    define("OneLoopDigits", d["digits"])
    a1, m2, devs = set(), set(), []
    for row in d["scan"]:
        f1, e1 = rational_times(float(mp.mpf(row["A1_over_pi"])), row["Nc"])
        f2, e2 = rational_times(float(mp.mpf(row["M2_limit_over_pi"])), row["Nc"])
        a1.add(f1); m2.add(f2)
        devs += [float(abs(mp.mpf(row["A1_over_pi"]) / (f1 * row["Nc"]) - 1)),
                 float(abs(mp.mpf(row["M2_limit_over_pi"]) / (f2 * row["Nc"]) - 1))]
    assert len(a1) == 1 and len(m2) == 1, "one-loop coefficients are not universal multiples of pi N_c"
    a1, m2 = a1.pop(), m2.pop()
    define("OneLoopAone", tex_fraction(a1, "\\pi\\Nc"))
    define("OneLoopMtwo", tex_fraction(m2, "\\pi\\Nc"))
    worst = max(devs)
    define("OneLoopAgreementDigits", d["stored_digits"] if worst == 0 else min(d["stored_digits"], int(math.floor(-math.log10(worst)))))
    define("OneLoopScanNc", listing(sorted({str(r["Nc"]) for r in d["scan"]})))
    define("OneLoopScanNf", listing(sorted({str(r["nf"]) for r in d["scan"]})))
    bench = next(a for a in d["angles"] if abs(float(mp.mpf(a["cos"])) - float(benchmark.COS_THETA)) < 1e-5)
    define("OneLoopAzero", f"{float(mp.mpf(bench['A0_re'])):.3f}{float(mp.mpf(bench['A0_im'])):+.3f}\\,i")
    # V1 at threshold: M1^(0) -> 2 i m_A^2, so 2 Re(M1^(1)/M1^(0)) -> (pi/2) A1 + Im A0
    v1_threshold = math.pi / 2 * float(mp.mpf(bench["A1"])) + float(mp.mpf(bench["A0_im"]))
    v = list(csv.DictReader(open(RESULTS / "virtual_correction.csv")))
    assert abs(float(v[0]["V1"]) - v1_threshold) < 0.05, "V1 at the smallest beta does not approach its threshold value"
    define("VoneThreshold", f"{v1_threshold:.2f}")
    return a1, m2


# ---------------------------------------------------------------- leading logarithms
def leading(a1, m2):
    d = json.load(open(RESULTS / "leading_logs.json"))
    nc, nf, lmu, ln2, pi = sp.symbols("nc n_f LMU LN2 PI")
    coeff = {k: sp.sympify(eval(v, {**sp.__dict__})) for k, v in d["exact_units_pi_mA2"].items()}
    names = {s.name: s for e in coeff.values() for s in e.free_symbols}
    Nc, Nf = names.get("nc", nc), names.get("n_f", nf)
    beta0 = (11 * Nc - 2 * Nf) / sp.Integer(3)

    def poly_tex(expr):
        """Nc^a (...) forms for the top coefficients, checked by expansion."""
        expr = sp.factor(expr)
        return expr

    c13, c12, c22, c21 = coeff["1,3"], coeff["1,2"], coeff["2,2"], coeff["2,1"]
    # closed forms in the text, each verified against the reduction
    forms = {"LeadOneThree": (c13, sp.Rational(-28, 3) * Nc**2, r"-\tfrac{28}{3}\pi\Nc^2"),
             "LeadOneTwo": (c12, sp.Rational(-14, 3) * Nc * (11 * Nc - 2 * Nf), r"-\tfrac{14}{3}\pi\Nc(11\Nc-2\nf)"),
             "LeadTwoTwo": (c22, -2 * Nc**2, r"-2\pi\Nc^2"),
             "LeadTwoOne": (c21, sp.Rational(-4, 3) * Nc * (4 * Nc - Nf), r"-\tfrac{4}{3}\pi\Nc(4\Nc-\nf)")}
    for name, (expr, form, tex) in forms.items():
        assert sp.expand(expr - form) == 0, name
        define(name, tex)
    # relation to the one-loop logarithm: c13 L^3 + c12 L^2 = a1 Nc L (x L^2 + y L)
    x = sp.simplify(c13 / (a1 * Nc)); y = sp.simplify(c12 / (a1 * Nc))
    assert sp.simplify(y + 2 * beta0) == 0, "the L^2 coefficient is no longer -2 beta0 times the one-loop log"
    define("LogRelationCubic", sp.latex(x).replace("nc", r"\Nc").replace(r"\frac", r"\tfrac"))
    define("MtwoRelation", sp.latex(sp.simplify(c22 / (m2 * Nc))).replace("nc", r"\Nc"))
    # the calL^1 coefficient of M1: (pi Nc / 9)[Nc(...) - nf(...)]
    c11 = sp.expand(coeff["1,1"] * 9 / Nc)
    mono = [(1, ""), (names.get("LN2", ln2), r"\ln2"), (names.get("LN2", ln2)**2, r"\ln^22"),
            (names.get("PI", pi)**2, r"\pi^2"), (sp.I * names.get("PI", pi), r"\,i\pi"), (names.get("LMU", lmu), r"L_\mu")]

    def bracket(poly):
        poly = sp.expand(poly)
        out, rest = "", poly
        for m, tex in mono:
            c = poly.coeff(m) if m != 1 else poly.subs({s: 0 for s in poly.free_symbols}).subs(sp.I, 0)
            if m == 1:
                c = sp.expand(poly).as_independent(*poly.free_symbols, as_Add=True)[0]
                c = sp.re(c)
            elif m == sp.I * names.get("PI", pi):
                c = sp.expand(poly).coeff(names.get("PI", pi)).coeff(sp.I)
            c = sp.nsimplify(c)
            if c == 0:
                continue
            out += ("+" if c > 0 and out else ("-" if c < 0 else "")) + (str(abs(c)) if (abs(c) != 1 or not tex) else "") + tex
            rest -= c * m
        assert sp.expand(rest) == 0, f"unexpected terms in the calL^1 coefficient: {sp.expand(rest)}"
        return out
    nc_part = c11.coeff(Nc); nf_part = -c11.coeff(Nf)
    assert sp.expand(c11 - Nc * nc_part + Nf * nf_part) == 0
    define("LeadOneOnePrefactor", r"\frac{\pi\Nc}{9}")
    define("LeadOneOneNc", bracket(nc_part))
    define("LeadOneOneNf", bracket(nf_part))
    define("LeadCheckedCount", d["coefficients_checked"])
    assert d["angle_only_through_calL_and_even"]


# ---------------------------------------------------------------- endpoint constants
def constants():
    d = json.load(open(RESULTS / "mb_constants.json"))
    define("ConstantCount", d["count"])
    define("ConstantCountWord", WORDS[d["count"]])
    define("PSLQDigits", d["digits"])
    define("PSLQResidual", sci(d["max_residual"], 1))
    rows = list(csv.DictReader(open(DATA / "constants" / "mb_constants.csv")))
    basis = set()
    for r in rows:
        expr = sp.expand(sp.sympify(r["reduction"], locals={"polylog": sp.polylog, "zeta": sp.zeta}))
        for term in sp.Add.make_args(expr):
            basis.add(term.as_coeff_Mul()[1])
    define("BasisCount", len(basis))
    define("BasisCountWord", WORDS[len(basis)])


# ---------------------------------------------------------------- normalization checks
def normalization():
    d = json.load(open(RESULTS / "normalization.json"))
    ratios = [float(mp.mpf(r)) for r in d["bornsqr_ratio"]]
    assert max(ratios) - min(ratios) < 1e-20
    define("BornSqrRatio", frac(ratios[0]))
    define("FeynmanDigits", int(math.floor(-math.log10(max(d["feynman_deviation"])))))
    imag = [float(mp.mpf(r["im"])) for r in d["nf_over_born"]]
    define("CataniImagDigits", int(math.floor(-math.log10(max(abs(v / (-2 * math.pi / 3) - 1) for v in imag) + 1e-300))))
    # one-loop files: dependence on m_A and mu_R only through their ratio after adding the scale term
    raw = [float(mp.mpf(r["raw_dev"])) for r in d["scaling"]]
    fixed = [float(mp.mpf(r["fixed_dev"])) for r in d["scaling"]]
    define("ScalingRawDev", sci(max(raw), 1))
    define("ScalingFixedDev", sci(max(fixed), 1))
    # two-loop scale dependence
    a = [float(mp.mpf(r["a"])) for r in d["rg"]]
    a_int = round(sum(a) / len(a))
    assert a_int == 3
    define("RGCoefficient", a_int)
    define("RGCoefficientPublished", a_int * FILE_FACTOR)
    define("RGaDev", sci(max(abs(v / a_int - 1) for v in a), 1))
    bdev = [abs(float(mp.mpf(r["b"])) / float(mp.mpf(r["b_expected"])) - 1) for r in d["rg"]]
    define("RGbDev", sci(max(bdev), 1))
    second = [abs(float(mp.mpf(r["second"])) - 1) for r in d["rg"]]
    define("RGSecondDev", sci(max(second), 1))
    define("RGSets", listing(sorted({f"$({r['Nc']},{r['nf']},{r['mA']}\\,\\mathrm{{GeV}})$" for r in d["rg"]})))


# ---------------------------------------------------------------- kinematics, cross section, virtual corrections
def physics():
    for b, tag in ((float(RECOMMENDED_BETA), "PointTwo"), (float(EDGE_BETA), "PointThree")):
        define(f"ThresholdPercent{tag}", f"{100 * (1 / math.sqrt(1 - b * b) - 1):.2f}")
        define(f"SqrtS{tag}", f"{2 * benchmark.M_A / math.sqrt(1 - b * b):.1f}")
    lo = list(csv.DictReader(open(RESULTS / "lo_cross_section.csv")))
    row = next(r for r in lo if abs(float(r["beta"]) - 0.1) < 1e-12)
    assert float(row["alpha_s"]) == benchmark.ALPHA_S and float(row["mA_GeV"]) == benchmark.M_A
    define("SigmaBeta", "0.1")
    define("SigmaLO", sci(float(row["sigma_hat_pb"])))
    v = list(csv.DictReader(open(RESULTS / "virtual_correction.csv")))
    a_s = benchmark.ALPHA_S / (4 * math.pi)
    v1 = [a_s * float(r["V1"]) for r in v]; v2 = [a_s**2 * float(r["V2"]) for r in v]
    define("VoneMin", f"\\ensuremath{{{min(v1):.2f}}}"); define("VoneMax", f"\\ensuremath{{{max(v1):.2f}}}")
    define("VtwoMin", f"\\ensuremath{{{min(v2):.2f}}}"); define("VtwoMax", f"\\ensuremath{{{max(v2):.2f}}}")
    define("VirtualBetaMin", f"\\ensuremath{{{power_of_ten(float(v[0]['beta']))}}}")
    define("VirtualBetaMax", f"{float(v[-1]['beta']):g}")
    # accuracy of the angular quadrature of scripts/virtual_correction.py, against 32 points
    from ggaa.crosssection import VirtualCorrections
    nodes = int(re.search(r'"--nodes", type=int, default=(\d+)', (ROOT / "scripts" / "virtual_correction.py").read_text()).group(1))
    vc, saved = VirtualCorrections(), mp.mp.dps
    mp.mp.dps = 30
    dev1 = dev2 = 1e-30
    betas = ("0.01", "0.1", "0.3")
    for b in betas:
        a1, a2 = vc.angle_averaged(mp.mpf(b), nodes=nodes)
        r1, r2 = vc.angle_averaged(mp.mpf(b), nodes=32)
        dev1 = max(dev1, float(abs(a1 / r1 - 1))); dev2 = max(dev2, float(abs(a2 / r2 - 1)))
    mp.mp.dps = saved
    words = {4: "four", 6: "six", 8: "eight", 10: "ten", 12: "twelve", 16: "sixteen"}
    define("QuadNodes", words.get(nodes, str(nodes))); define("QuadNodesRef", "32")
    define("QuadBetas", "$" + "$, $".join(betas[:-1]) + "$ and $" + betas[-1] + "$")
    define("QuadVoneDev", f"\\ensuremath{{10^{{{math.ceil(math.log10(dev1))}}}}}")
    define("QuadVtwoDev", f"\\ensuremath{{10^{{{math.ceil(math.log10(dev2))}}}}}")


# ---------------------------------------------------------------- timing and sizes
def timing():
    rows = list(csv.DictReader(open(DATA / "timing" / "timing_benchmark.csv")))
    first = sorted(float(r["FirstCallSpeedup"]) for r in rows)
    warm = sorted(float(r["WarmSpeedup"]) for r in rows)
    median = lambda xs: xs[len(xs) // 2]
    define("TimingFirstMin", f"{first[0]:.2f}"); define("TimingFirstMax", f"{first[-1]:.2f}")
    define("TimingFirstMedian", f"{median(first):.2f}")
    define("TimingWarmMin", f"{warm[0]:.2f}"); define("TimingWarmMax", f"{warm[-1]:.2f}")
    define("TimingWarmMedian", f"{median(warm):.2f}")
    define("TimingPrecision", rows[0]["WorkingPrecision"])
    define("TimingRepeats", WORDS[int(rows[0]["WarmRepeats"])])
    define("TimingCos", rows[0]["CosTheta"])
    define("TimingBetas", listing([f"${r['Beta']}$" for r in rows]))
    m2, mu2, nc = (int(x) for x in rows[0]["Parameters"].strip("{}").split(","))
    define("TimingMA", int(round(math.sqrt(m2)))); define("TimingMuR", int(round(math.sqrt(mu2)))); define("TimingNc", nc)
    run = json.load(open(DATA / "timing" / "timing_run.json"))["representation_bytes"]
    define("TimingExactMiB", f"{run['unexpanded'] / 2**20:.2f}")
    define("TimingSeriesMiB", f"{run['beta6'] / 2**20:.2f}")
    define("TimingRepresentationRatio", f"{run['unexpanded'] / run['beta6']:.1f}")
    published = [INPUTS / "Finite_2Loop_P1", INPUTS / "Finite_2Loop_P2"]
    if not all(p.exists() for p in published):
        sys.exit("make inputs first: the size of the published two-loop files is read from inputs/")
    published_bytes = sum(p.stat().st_size for p in published)
    series_bytes = sum(p.stat().st_size for p in SERIES.glob("M*_nf*.m"))
    define("PublishedTwoLoopMB", f"{published_bytes / 1e6:.0f}")
    define("SeriesKB", f"{series_bytes / 1e3:.0f}")
    define("SizeRatio", f"{published_bytes / series_bytes:.0f}")


def loop_sizes():
    """|a_s^n M_i^(n)| / |M_i^(0)| at beta = 1/10 and the benchmark angle."""
    from ggaa.born import born
    from ggaa.oneloop import OneLoop
    from ggaa.series import ThresholdSeries
    if not (INPUTS / "Finite_1Loop_P1").exists():
        sys.exit("make inputs first")
    mp.mp.dps = 60
    a_s = benchmark.ALPHA_S / (4 * math.pi)
    b = mp.mpf(1) / 10
    c = mp.mpf(benchmark.COS_THETA.numerator) / benchmark.COS_THETA.denominator
    tree = abs(born(1, b, c, benchmark.M_A))
    define("LoopSizeBeta", "0.1")
    define("MoneOneLoopRatio", f"{float(a_s * abs(OneLoop()(1, b, c)) / tree):.1f}")
    define("MoneTwoLoopRatio", f"{float(a_s**2 * abs(ThresholdSeries()(1, b, c)) / tree):.1f}")
    define("MtwoOneLoopRatio", f"{float(a_s * abs(OneLoop()(2, b, c)) / abs(born(2, b, c, benchmark.M_A))):.1f}")


def structure():
    """Double parton scattering check, logarithms in V2, top-mass factor, quark-loop term."""
    import json
    import sympy as sp
    define("BenchMt", f"{benchmark.M_TOP:g}")
    dps = json.load(open(RESULTS / "dps_check.json"))
    betas = sorted({p["beta"] for p in dps["points"]})
    define("DPSMaxDeviation", sci(dps["max_relative_difference"], 1))
    # one-loop gg -> HH through two CP-even vertices (scripts/check_dps_hh.py, data/dps_hh)
    hh = json.load(open(RESULTS / "dps_hh.json"))
    define("HHDpsMaxDev", sci(hh["max_dev_HH"], 1))
    define("HHMtwoMax", sci(hh["max_M2_over_M1"], 1))
    define("HHBetaMin", f"{hh['beta_min']:g}")
    define("HHBetaMax", f"{hh['beta_max']:g}")
    define("HHReproduceDigits", str(hh["aa_reproduction_digits"]))
    # the 2019 renormalisation and the June 2019 pole record (scripts/check_renorm_2019.py)
    rn = json.load(open(RESULTS / "renorm_2019.json"))
    define("RenormPointsWord", {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten"}[rn["pole_record_n_points"]])
    define("DPSBetaMin", f"{betas[0]:g}")
    define("DPSBetaMax", f"{betas[-1]:g}")
    define("DPSNcValues", listing([str(n) for n in sorted({p["Nc"] for p in dps["points"]})]))
    define("DPSNfValues", listing([str(n) for n in sorted({p["nf"] for p in dps["points"]})]))
    logs = json.load(open(RESULTS / "rate_logs.json"))
    define("VtwoLogFit", f"{logs['fit_c']:.6f}")
    define("VtwoLogExact", f"{logs['expected_c']:.6f}")
    define("VtwoLogRelDiff", sci(logs["relative_difference"], 1))
    a_s = benchmark.ALPHA_S / (4 * math.pi)
    define("VtwoLogScaled", f"\\ensuremath{{{a_s**2 * logs['expected_c']:.4f}}}")
    define("OneLoopScaleShift", f"{2 * (11 * benchmark.N_C - 2 * benchmark.N_F) / 3 * math.log(benchmark.M_A**2):.1f}")
    exps = sorted(int(b.split("e-")[1]) for b in logs["betas"])
    define("VtwoBetaSmall", f"10^{{-{exps[-1]}}}")
    define("VtwoBetaLarge", f"10^{{-{exps[0]}}}")
    top = json.load(open(RESULTS / "top_mass.json"))
    masses = sorted(top["threshold_ratio"], key=int)
    for name, mA in zip(("A", "B", "C"), masses):
        define(f"TopMassMA{name}", mA)
        define(f"TopMassRatio{name}", f"{top['threshold_ratio'][mA]:.3f}")
    worst = max(top["threshold_ratio"].values())
    define("TopMassPercentMax", f"{(worst - 1) * 100:.1f}")
    for name, mA in zip(("A", "B", "C"), masses):
        define(f"TopMassDPS{name}", f"{top['dps_vertex_ratio'][mA]:.3f}")
    # size of the logarithmic terms of V2 at fixed angle (scripts/check_log_share.py)
    share = json.load(open(RESULTS / "log_share.json"))
    small, large = sorted(share["points"], key=lambda r: r["beta"])
    assert abs(share["cos"] - 1 / 3) < 1e-12
    define("LogShareBetaSmall", f"{small['beta']:g}")
    define("LogShareBetaLarge", f"{large['beta']:g}")
    define("LogShareSmall", f"{100 * small['log_share']:.0f}")
    define("LogShareLarge", f"{100 * large['log_share']:.0f}")
    define("SquareShareSmall", f"{100 * abs(small['square_share']):.0f}")
    define("SquareShareLarge", f"{100 * abs(large['square_share']):.0f}")
    # leading-power structure of V2 (scripts/check_rate_logs.py) and the power corrections at beta = 0.1
    rl = json.load(open(RESULTS / "rate_logs.json"))
    define("VtwoAngleSpread", sci(rl["angle_spread_v0_v1"], 1))
    define("VoneScaleDev", sci(max(rl["v1_scale_rel_diff"], 1e-30), 1))
    define("VoneFormulaDev", sci(max(rl["v1_formula_rel_diff"], 1e-30), 1))
    define("VzeroBench", f"{rl['fit_a']:.0f}")
    define("VzeroTwoMu", f"{rl['scale_fit_2muR']['v0']:.0f}")
    define("VoneFit", f"{rl['fit_b']:.1f}")
    define("VtwoPeakEll", f"{rl['ell_peak']:.2f}")
    define("VtwoPeakBeta", f"{rl['beta_peak']:.4f}")
    lead_dev = abs(large["V2"] - large["leading_power"]) / abs(large["V2"])
    define("LeadPowerDevLarge", f"{100 * lead_dev:.1f}")
    # velocity below which the two-loop term of M_1 exceeds the one-loop term (check_amplitude_crossover.py)
    cross = json.load(open(RESULTS / "amplitude_crossover.json"))
    define("CrossoverBeta", sci(cross["beta_crossover"], 2))
    # comparison with the 2019 archive (scripts/check_archive_2019.py)
    arch = json.load(open(RESULTS / "archive_2019.json"))
    define("ArchiveOneLoopDev", sci(float(mp.mpf(arch["one_loop"]["max_rel_dev_from_corrected"])), 1))
    define("ArchiveTwoLoopDev", sci(float(mp.mpf(arch["two_loop"]["real_part_rel_dev"])), 1))
    res = float(mp.mpf(arch["two_loop"]["residual_over_32CF_G0"]))
    assert abs(res - 1) < 1e-15
    a_s = benchmark.ALPHA_S / (4 * math.pi)
    oj_shift = a_s**2 * 2 * float(mp.mpf(arch["two_loop"]["residual_over_M0"])) / FILE_FACTOR
    define("ArchiveOJShift", f"{oj_shift:.2f}")
    # the published one-loop O_G-O_J remainders, Eqs. (3.6)-(3.7) of arXiv:1909.08993, over the Born
    x, c, nf = sp.symbols("x c n_f", positive=True)
    b = (1 - x) / (1 + x)
    z = (1 + b**2 + 2 * b * c) / (1 - b**2)
    born1 = sp.I * (1 + x) ** 2 / (2 * x)
    born2 = sp.I * (x - z) * (x**2 + 1) * (x * z - 1) / (4 * x * z * (x**2 - x * z + 1))
    gj1 = -24 * sp.I * nf * (1 + x) ** 2 / x
    gj2 = -12 * sp.I * nf * (1 - x**2) ** 2 * (1 + x**2) * (1 - c**2) / (x * ((1 + x**2) ** 2 - (1 - x**2) ** 2 * c**2))
    r1, r2 = sp.simplify(gj1 / born1 / nf), sp.simplify(gj2 / born2 / nf)
    assert r1 == r2 and r1.is_Integer, (r1, r2)
    define("GJRatio", str(r1))


# ---------------------------------------------------------------- forward limit of the two-loop |J_z| = 2 amplitude
def jz2():
    """Macros for Sec. IV D from results/jz2_forward_points.jsonl (scripts/check_jz2_forward.py).

    For each colour and flavour set and velocity, the two-loop M_2 at 1 - c = 10^-5, 10^-6, 10^-7
    (140 digits) is fitted exactly by A u^2 + B u + C, u = ln(1 - c); a2 = 4 A is the coefficient of
    l^2, l = ln(2 p_T / m).  It is compared with -4 N_c K1, K1 the forward one-loop M_2 / m^2.
    """
    from ggaa.oneloop import OneLoop
    rows = [json.loads(l) for l in (RESULTS / "jz2_forward_points.jsonl").read_text().splitlines() if l.strip()]
    exps, dps, check = (5, 6, 7), 140, 180
    sets = sorted({(r["Nc"], r["nf"]) for r in rows})
    saved, amp = mp.mp.dps, OneLoop()
    worst, digits, betas, k1 = 0.0, 1e9, set(), {}
    for Nc, nf in sets:
        for b in sorted({r["beta"] for r in rows if (r["Nc"], r["nf"]) == (Nc, nf)}, key=float):
            pick = {r["e"]: r for r in rows if (r["Nc"], r["nf"], r["beta"], r["dps"]) == (Nc, nf, b, dps)}
            if not all(e in pick for e in exps):
                continue
            mp.mp.dps = 60
            us = [-e * mp.log(10) for e in exps]
            vs = [mp.mpc(mp.mpf(pick[e]["re"]), mp.mpf(pick[e]["im"])) for e in exps]
            A = mp.lu_solve(mp.matrix([[u**2, u, 1] for u in us]), mp.matrix(vs))[0]
            mp.mp.dps = 120
            K1 = (amp(2, mp.mpf(b), 1 - mp.mpf(10)**-18, mA=125, muR=125, Nc=Nc, nf=nf) / 125**2).real
            worst = max(worst, float(abs((4 * A).real / (-4 * Nc * K1) - 1)))
            k1.setdefault((Nc, nf), []).append(float(K1))
            betas.add(b)
            hi = [r for r in rows if (r["Nc"], r["nf"], r["beta"], r["e"], r["dps"]) == (Nc, nf, b, max(exps), check)]
            if hi:
                mp.mp.dps = 60
                v1 = mp.mpc(mp.mpf(hi[0]["re"]), mp.mpf(hi[0]["im"])); v0 = vs[-1]
                digits = min(digits, float(-mp.log10(abs(v1 - v0) / abs(v1))))
    mp.mp.dps = saved
    exponent = math.floor(math.log10(worst)); mantissa = math.ceil(worst / 10**exponent - 1e-9)
    if mantissa == 10:
        mantissa, exponent = 1, exponent + 1
    define("JzMaxDev", f"\\ensuremath{{{mantissa}\\times10^{{{exponent}}}}}")
    define("JzSets", " and ".join(f"$({a},{b})$" for a, b in sets))
    define("JzBetaMin", f"{min(betas, key=float)}"); define("JzBetaMax", f"{max(betas, key=float)}")
    define("JzDigits", str(dps)); define("JzCheckDigits", str(check))
    define("JzCheckAgreement", str(int(math.floor(digits))))
    define("JzExponents", "$10^{-5}$, $10^{-6}$ and $10^{-7}$")
    ratio = max(max(v) / min(v) for v in k1.values())
    define("JzKOneRatio", f"{ratio:.1f}")


def main():
    parameters()
    loop_sizes()
    err = validation()
    forward(err)
    exact_check()
    a1, m2 = one_loop()
    leading(a1, m2)
    constants()
    normalization()
    physics()
    timing()
    structure()
    jz2()
    target = ROOT / "manuscript" / "numbers.tex"
    lines = ["% Generated by scripts/make_numbers.py from the data and results folders. Do not edit."]
    lines += [f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in sorted(macros.items())]
    target.write_text("\n".join(lines) + "\n")
    print(f"wrote manuscript/numbers.tex ({len(macros)} macros)")


if __name__ == "__main__":
    main()
