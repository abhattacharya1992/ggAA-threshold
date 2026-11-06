"""Endpoint constants of the threshold expansion.

The construction leaves thirteen constants

    MB[p, q, a, w] = a/(p-1)! * Integral_0^1 dt (-ln t)^(p-1) / (1 - a t) * phi_q(w t),

with phi_1(z) = -ln(1 - z) and phi_q(z) = Li_q(z) for q >= 2.  All of them reduce to
the standard basis of alternating sums of weight three and four.  The reductions below
were found with PSLQ at 110 digits; `identify` repeats that search.
"""

from fractions import Fraction

import mpmath as mp
import sympy as sp

ZETA3 = sp.zeta(3)
LN2 = sp.log(2)
LI4_HALF = sp.polylog(4, sp.Rational(1, 2))
PI = sp.pi

HALF = Fraction(1, 2)

REDUCTIONS = {
    (1, 2, -1, 1): sp.Rational(5, 8) * ZETA3 - PI**2 * LN2 / 6,
    (1, 2, HALF, HALF): -ZETA3 / 4 + PI**2 * LN2 / 12 - LN2**3 / 6,
    (1, 2, HALF, 1): -ZETA3 + PI**2 * LN2 / 4,
    (2, 1, -1, -1): ZETA3 / 8,
    (2, 1, -1, 1): -sp.Rational(13, 8) * ZETA3 + PI**2 * LN2 / 4,
    (2, 1, 1, -1): ZETA3 - PI**2 * LN2 / 4,
    (2, 1, 1, HALF): -sp.Rational(5, 8) * ZETA3 + PI**2 * LN2 / 6,
    (2, 1, HALF, HALF): ZETA3 / 8 - LN2**3 / 6,
    (2, 1, HALF, 1): sp.Rational(13, 8) * ZETA3 - PI**2 * LN2 / 4,
    (2, 2, 1, -1): (-sp.Rational(71, 1440) * PI**4 - PI**2 * LN2**2 / 6 + LN2**4 / 6
                    + 4 * LI4_HALF + sp.Rational(7, 2) * ZETA3 * LN2),
    (2, 2, 1, HALF): (sp.Rational(19, 1440) * PI**4 - PI**2 * LN2**2 / 24 - LN2**4 / 24
                      - LI4_HALF - ZETA3 * LN2 / 4),
    (2, 2, HALF, HALF): (PI**4 / 1440 - PI**2 * LN2**2 / 24 + LN2**4 / 24
                         + ZETA3 * LN2 / 4),
    (2, 2, HALF, 1): (sp.Rational(7, 288) * PI**4 - PI**2 * LN2**2 / 8 - LN2**4 / 8
                      - 3 * LI4_HALF),
}


def _mp(x):
    x = Fraction(x)
    return mp.mpf(x.numerator) / x.denominator


def mb_integral(p, q, a, w):
    """MB[p, q, a, w] by numerical quadrature at the current mpmath precision."""
    a, w = _mp(a), _mp(w)
    phi = (lambda z: -mp.log(1 - z)) if q == 1 else (lambda z: mp.polylog(q, z))

    def integrand(t):
        return a * (-mp.log(t)) ** (p - 1) / (1 - a * t) * phi(w * t) / mp.factorial(p - 1)

    return mp.quad(integrand, [0, 0.25, 0.5, 0.75, 1])


def basis(weight):
    """Names and values of the product basis of the given weight (3 or 4)."""
    ln2, z3, li4 = mp.log(2), mp.zeta(3), mp.polylog(4, mp.mpf(1) / 2)
    if weight == 3:
        return (["zeta3", "pi^2 ln2", "ln2^3"],
                [z3, mp.pi**2 * ln2, ln2**3])
    if weight == 4:
        return (["pi^4", "pi^2 ln2^2", "ln2^4", "Li4(1/2)", "zeta3 ln2"],
                [mp.pi**4, mp.pi**2 * ln2**2, ln2**4, li4, z3 * ln2])
    raise ValueError("weight must be 3 or 4")


def identify(p, q, a, w):
    """PSLQ relation of MB[p, q, a, w] with the basis of weight p + q.

    Returns (value, {basis name: rational coefficient}, residual)."""
    value = mb_integral(p, q, a, w)
    names, values = basis(p + q)
    relation = mp.pslq([value] + values, maxcoeff=10**7, maxsteps=10**6)
    if relation is None:
        raise RuntimeError(f"no relation found for MB[{p},{q},{a},{w}]")
    coefficients = {n: Fraction(-c, relation[0]) for n, c in zip(names, relation[1:]) if c}
    residual = abs(relation[0] * value + sum(c * v for c, v in zip(relation[1:], values)))
    return value, coefficients, residual
