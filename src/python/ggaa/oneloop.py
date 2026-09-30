"""Direct evaluation of the published one-loop Type-II finite remainders.

The ancillary files Finite_1Loop_P1 and Finite_1Loop_P2 of arXiv:1909.08993 are
Wolfram Language expressions in x, z, Y = 1 + y, Z = 1 + z, M2 = m_A^2, muR2, N, Nf.
They are translated once into Python and evaluated with mpmath on the threshold map

    x = (1 - b)/(1 + b),  y = (1 + b^2 - 2 b c)/(1 - b^2),  z = (1 + b^2 + 2 b c)/(1 - b^2).

Near threshold the expressions cancel strongly; use a working precision well above
the number of digits wanted (see scripts/check_one_loop.py).

Scale dependence.  The published one-loop files do not depend on mu_R: their eight
Log[muR2] terms cancel, and the remainders change under a common rescaling of m_A and
mu_R, i.e. they contain ln m_A^2 where ln(m_A^2/mu_R^2) is required.  Renormalization-group
invariance (C_G proportional to a_s, Z_GG = Z_as, Catani's I operators) requires
dM^(1),fin/dln mu_R^2 = 2 beta0 M^(0), and the published two-loop remainders are consistent
with the one-loop remainders only if

    M^(1),fin(mu_R) = [published file] + 2 beta0 ln(mu_R^2 / GeV^2) M^(0),

i.e. the published files equal the remainders at mu_R = 1 GeV (scripts/check_normalization.py,
checks 4 and 5).  __call__ applies this term by default; scale_fix=False returns the file.
"""

import re

import mpmath as mp

from .born import born
from .paths import INPUTS

FILES = {1: "Finite_1Loop_P1", 2: "Finite_1Loop_P2"}


def _to_python(text):
    text = text.split("=", 1)[1] if text.lstrip().startswith("M_") else text
    text = text.strip().rstrip(";")
    text = text.replace("[", "(").replace("]", ")").replace("^", "**")
    text = re.sub(r"\bLog\(", "log(", text)
    text = re.sub(r"\bPolyLog\(", "polylog(", text)
    text = re.sub(r"\bPi\b", "pi", text)
    text = re.sub(r"\bI\b", "j", text)
    text = re.sub(r"(?<![\w.])(\d+)(?![\w.])", r"mpf(\1)", text)  # exact integers
    return "(" + text + ")"


class OneLoop:
    def __init__(self, directory=INPUTS):
        self._code = {i: compile(_to_python((directory / name).read_text()), name, "eval")
                      for i, name in FILES.items()}

    def __call__(self, i, beta, costh, mA=125, muR=125, Nc=3, nf=5, scale_fix=True):
        """M_i^{II,(1),fin} at the given point (couplings stripped), masses in GeV.

        With scale_fix=True (default) the term 2 beta0 ln(mu_R^2/GeV^2) M^(0) missing from the
        published files is added; scale_fix=False returns the published expression."""
        b, c = mp.mpf(beta), mp.mpf(costh)
        x = (1 - b) / (1 + b)
        y = (1 + b**2 - 2 * b * c) / (1 - b**2)
        z = (1 + b**2 + 2 * b * c) / (1 - b**2)
        names = dict(mpf=mp.mpf, log=mp.log, pi=mp.pi, j=mp.mpc(0, 1),
                     polylog=lambda n, a: mp.polylog(int(n), a),
                     x=x, z=z, Y=1 + y, Z=1 + z, M2=mp.mpf(mA) ** 2, muR2=mp.mpf(muR) ** 2,
                     N=mp.mpf(Nc), Nf=mp.mpf(nf))
        value = eval(self._code[i], names)
        if scale_fix:
            beta0 = (11 * mp.mpf(Nc) - 2 * mp.mpf(nf)) / 3
            value += 2 * beta0 * mp.log(mp.mpf(muR) ** 2) * born(i, b, c, mA)
        return value
