"""Read the Wolfram Language coefficient files of the series construction and bring
them to a minimal set of functions and constants.

After `reduce` every coefficient is a polynomial in the symbols below, with
coefficients that are rational functions of cos(theta), m_A^2 and N_c:

    L    = ln(beta_A)
    LM   = ln(1 - cos theta),  LP = ln(1 + cos theta)
    LI2  = Li2((1 + cos theta)/2)
    LMU  = ln(m_A^2 / mu_R^2)
    LN2  = ln 2,  PI = pi,  Z3 = zeta(3),  LI4 = Li4(1/2)

The endpoint constants MB[p,q,a,w] are replaced by their reductions in
`ggaa.constants`.  Dilogarithms of other angular arguments are mapped onto LI2 with
the reflection and Landen identities, valid for -1 < cos theta < 1.  With these
substitutions an exact cancellation becomes a polynomial identity.
"""

import re
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.parsing.mathematica import parse_mathematica
from sympy.printing.mathematica import mathematica_code

from . import constants

costh, m2, mu2, nc = sp.symbols("costh m2 mu2 nc")
L, LM, LP, LI2, LMU, LN2, PI, Z3, LI4 = sp.symbols("L LM LP LI2 LMU LN2 PI Z3 LI4")
TRANSCENDENTALS = (L, LM, LP, LI2, LMU, LN2, PI, Z3, LI4)

_LOG_M2, _LOG_MU2, _BETA = sp.symbols("_log_m2 _log_mu2 bta")

# One term of the n_f^2 coefficient of M_2 at beta^6 is stored as the logarithm of a
# complex power.  Its value is the principal branch, which is the branch used in the
# validation and the one that agrees with the unexpanded amplitude.
_PRINCIPAL_LOG = ("Log[4^(24*(7 + (6*I)*Pi))]",
                  "(168*Log[4] + I*(144*Pi*Log[4] - 200*Pi))")


def read(path):
    """Parse one coefficient file into a sympy expression."""
    text = Path(path).read_text().strip()
    text = re.sub(r"\(\*.*?\*\)", "", text, flags=re.S)
    text = text.replace(*_PRINCIPAL_LOG)
    # `beta` would be read as sympy's Beta function.
    text = re.sub(r"\bbeta\b", "bta", text)
    return parse_mathematica(text).subs(sp.Symbol("Lb"), L)


def _li2(z):
    """Li2 of the angular arguments that occur, in terms of LI2 = Li2((1+c)/2)."""
    li2_minus = PI**2 / 6 - (LM - LN2) * (LP - LN2) - LI2          # Li2((1-c)/2)
    known = [((1 + costh) / 2, LI2),
             ((1 - costh) / 2, li2_minus),
             ((costh - 1) / (costh + 1), -li2_minus - (LP - LN2) ** 2 / 2),
             ((costh + 1) / (costh - 1), -LI2 - (LM - LN2) ** 2 / 2)]
    for argument, value in known:
        if sp.cancel(z - argument) == 0:
            return value
    if z == sp.Rational(1, 2):
        return PI**2 / 12 - LN2**2 / 2
    raise ValueError(f"unexpected dilogarithm argument {z}")


def _polylog(n, z):
    n = int(n)
    if n == 2:
        return _li2(z)
    if n == 3 and z == sp.Rational(1, 2):
        return sp.Rational(7, 8) * Z3 - PI**2 * LN2 / 12 + LN2**3 / 6
    if n == 4 and z == sp.Rational(1, 2):
        return LI4
    raise ValueError(f"unexpected PolyLog[{n}, {z}]")


def _log(z):
    if z.is_Rational and z > 0:
        num, den = sp.factorint(z.p), sp.factorint(z.q)
        if set(num) <= {2} and set(den) <= {2}:
            return (num.get(2, 0) - den.get(2, 0)) * LN2
    for argument, value in ((1 - costh, LM), (1 + costh, LP)):
        if sp.cancel(z - argument) == 0:
            return value
    named = {m2: _LOG_M2, mu2: _LOG_MU2, m2 / mu2: LMU, _BETA: L}
    if z in named:
        return named[z]
    raise ValueError(f"unexpected logarithm argument {z}")


def _mb(p, q, a, w):
    key = (int(p), int(q), Fraction(str(a)), Fraction(str(w)))
    return constants.REDUCTIONS[key]


def reduce(expr):
    """Rewrite a parsed coefficient in the symbols of this module and expand it."""
    expr = expr.replace(sp.Function("MB"), _mb)
    expr = expr.replace(sp.Function("PolyLog"), _polylog).replace(sp.polylog, _polylog)
    zeta = lambda n: Z3 if n == 3 else sp.zeta(n)
    expr = expr.replace(sp.Function("Zeta"), zeta).replace(sp.zeta, zeta)
    expr = expr.replace(sp.log, _log).subs(sp.pi, PI)
    expr = sp.expand(expr.subs(_LOG_M2, _LOG_MU2 + LMU))
    if expr.has(_LOG_MU2):
        raise ValueError("ln(mu_R^2) does not enter only through ln(m_A^2/mu_R^2)")
    return expr


def collect(expr):
    """Map each transcendental monomial to its (cancelled) rational coefficient."""
    groups = defaultdict(int)
    for term in sp.Add.make_args(sp.expand(expr)):
        monomial, coefficient = sp.Integer(1), sp.Integer(1)
        for factor in sp.Mul.make_args(term):
            if factor.free_symbols & set(TRANSCENDENTALS):
                monomial *= factor
            else:
                coefficient *= factor
        groups[monomial] += coefficient
    groups = {m: sp.factor(sp.cancel(sp.together(v))) for m, v in groups.items()}
    return {m: v for m, v in groups.items() if v != 0}


def is_zero(expr):
    return sp.expand(sp.numer(sp.together(expr))) == 0


def calL_shift():
    """L = calL - shift, with calL = L + (LM + LP)/2 + LN2 + i pi/2 = (1/2) ln[-4 beta^2 (1-c^2)]."""
    return (LM + LP) / 2 + LN2 + sp.I * PI / 2


def to_wolfram(expr):
    """Wolfram Language input form, one monomial per line (each line but the last ends in +)."""
    placeholders = {s: sp.Symbol(f"X{s.name}X") for s in TRANSCENDENTALS}
    names = {"XLX": "Log[beta]", "XLMX": "Log[1 - costh]", "XLPX": "Log[1 + costh]",
             "XLI2X": "PolyLog[2, (1 + costh)/2]", "XLMUX": "Log[m2/mu2]",
             "XLN2X": "Log[2]", "XPIX": "Pi", "XZ3X": "Zeta[3]", "XLI4X": "PolyLog[4, 1/2]"}
    terms = collect(expr)
    if not terms:
        return "0"
    lines = []
    for monomial in sorted(terms, key=lambda m: (-sp.Poly(m, *TRANSCENDENTALS).total_degree(), str(m))):
        mono = mathematica_code(monomial.subs(placeholders))
        for key, value in names.items():
            mono = mono.replace(key, value)
        coefficient = mathematica_code(terms[monomial])
        lines.append(f"({coefficient})" if monomial == 1 else f"({coefficient})*{mono}")
    return " +\n".join(lines)
