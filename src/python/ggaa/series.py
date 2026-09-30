"""Evaluate the threshold expansion of the two-loop Type-II finite remainders.

The files data/series/M{i}_nf{k}.m hold the coefficient of n_f^k of the tensor
coefficient M_i^{II,(2),fin} as a list {c0, c2, c4, c6}:

    M_i,k = c0 + beta^2 c2 + beta^4 c4 + beta^6 c6 + O(beta^8),

where each c_n is a polynomial in ln(beta).  Odd powers and inverse powers of beta
vanish identically (scripts/reduce_series.py, `make series`).  The complete coefficient is
M_i = M_i,0 + n_f M_i,1 + n_f^2 M_i,2.  Normalization as in ggaa.normalization: the
amplitude is C_G^2 [M^(0) + a_s M^(1) + a_s^2 M^(2)] delta^{ab}, which is the expansion of the
published two-loop files of arXiv:1909.08993 divided by 64.
"""

import re
from pathlib import Path

import mpmath as mp
import sympy as sp
from sympy.parsing.mathematica import parse_mathematica

from .paths import SERIES

ORDERS = (0, 2, 4, 6)
_ARGS = sp.symbols("bta costh m2 mu2 nc")


def read_series_file(path):
    """Return the coefficients {c0, c2, c4, c6} of one series file as sympy expressions."""
    text = Path(path).read_text()
    text = re.sub(r"\(\*.*?\*\)", "", text, flags=re.S)
    text = re.sub(r"\bbeta\b", "bta", text)
    parsed = parse_mathematica(text.strip())
    parsed = parsed.replace(sp.Function("PolyLog"), sp.polylog)
    parsed = parsed.replace(sp.Function("Zeta"), sp.zeta)
    coefficients = tuple(parsed)
    if len(coefficients) != len(ORDERS):
        raise ValueError(f"{path}: expected {len(ORDERS)} coefficients")
    return coefficients


class ThresholdSeries:
    """The saved expansion of M_1 and M_2 through beta^6."""

    def __init__(self, directory=SERIES):
        self._functions = {}
        for i in (1, 2):
            for k in (0, 1, 2):
                coefficients = read_series_file(directory / f"M{i}_nf{k}.m")
                self._functions[i, k] = [sp.lambdify(_ARGS, c, modules="mpmath")
                                         for c in coefficients]

    def component(self, i, k, beta, costh, mA=125, muR=125, Nc=3, order=6):
        """Coefficient of n_f^k in M_i, truncated after beta^order (order = 4 or 6)."""
        if order not in (4, 6):
            raise ValueError("order must be 4 or 6")
        beta, costh = _mpf(beta), _mpf(costh)
        args = (beta, costh, _mpf(mA) ** 2, _mpf(muR) ** 2, _mpf(Nc))
        return sum(beta**n * f(*args)
                   for n, f in zip(ORDERS, self._functions[i, k]) if n <= order)

    def __call__(self, i, beta, costh, mA=125, muR=125, Nc=3, nf=5, order=6):
        """M_i^{II,(2),fin} through beta^order for n_f light flavours."""
        nf = _mpf(nf)
        return sum(nf**k * self.component(i, k, beta, costh, mA, muR, Nc, order)
                   for k in (0, 1, 2))


def _mpf(x):
    """Exact rationals given as strings such as '1/3' are converted without rounding."""
    if isinstance(x, str) and "/" in x:
        num, den = x.split("/")
        return mp.mpf(num) / mp.mpf(den)
    return mp.mpf(x)
