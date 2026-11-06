"""Physical quantities for the Type-II contribution to g g -> A A.

Helicity sums.  T_1 contributes only to equal gluon helicities and T_2 only to opposite
ones, each with a unit-modulus contraction, so

    sum_{hel, colour} |A|^2 = 2 (N_c^2 - 1) C_G^4 ( |M_1|^2 + |M_2|^2 ).

Partonic cross section.  Averaging over 2 x 2 helicities and (N_c^2 - 1)^2 colours, and
including 1/2 for identical particles,

    d sigma / d cos(theta) = beta C_G^4 ( |M_1|^2 + |M_2|^2 ) / ( 128 pi s (N_c^2 - 1) ),

integrated over -1 < cos(theta) < 1.  With O_G = eps G G as in arXiv:1909.08993,
C_G = -2^(5/4) a_s G_F^(1/2) cot(beta_H), a_s = alpha_s / (4 pi).

Virtual corrections.  With M_i = C_G^2 [M_i^(0) + a_s M_i^(1) + a_s^2 M_i^(2)],

    sum|M|^2 = C_G^4 sum_i |M_i^(0)|^2 [ 1 + a_s V1 + a_s^2 V2 ],
    V1 = 2 Re sum_i M_i^(0)* M_i^(1) / sum_i |M_i^(0)|^2,
    V2 = [ 2 Re sum_i M_i^(0)* M_i^(2) + sum_i |M_i^(1)|^2 ] / sum_i |M_i^(0)|^2,

with the finite remainders M^(1), M^(2) (Catani subtraction) in one normalization
(ggaa.normalization).  V2 does not include the O_G-O_J term with the physical C_J.
V1 and V2 are the finite virtual parts of an NLO and NNLO calculation.  They depend on the
infrared subtraction scheme and are not observables by themselves: real radiation and
mass factorization are needed for a physical cross section.
"""

import mpmath as mp

from .born import born
from .oneloop import OneLoop
from .series import ThresholdSeries

GEV2_TO_PB = mp.mpf("0.3893793721e9")    # 1 GeV^-2 in pb
G_FERMI = mp.mpf("1.1663787e-5")         # GeV^-2


def c_g(alpha_s, cot_beta=1, g_fermi=G_FERMI):
    """Wilson coefficient C_G in GeV^-1 (leading order in alpha_s)."""
    a_s = mp.mpf(alpha_s) / (4 * mp.pi)
    return -mp.mpf(2) ** (mp.mpf(5) / 4) * a_s * mp.sqrt(g_fermi) * mp.mpf(cot_beta)


def born_sum(beta, costh, mA=125):
    """|M_1^(0)|^2 + |M_2^(0)|^2 (GeV^4, couplings stripped)."""
    return sum(abs(born(i, beta, costh, mA)) ** 2 for i in (1, 2))


def dsigma_dcos(beta, costh, mA=125, alpha_s=0.118, cot_beta=1, Nc=3):
    """LO partonic d sigma / d cos(theta) of the Type-II amplitude, in pb."""
    beta = mp.mpmathify(beta)
    s = 4 * mp.mpf(mA) ** 2 / (1 - beta**2)
    value = beta * c_g(alpha_s, cot_beta) ** 4 * born_sum(beta, costh, mA) / (128 * mp.pi * s * (Nc**2 - 1))
    return value * GEV2_TO_PB


def sigma(beta, mA=125, alpha_s=0.118, cot_beta=1, Nc=3):
    """LO partonic cross section in pb."""
    return mp.quad(lambda c: dsigma_dcos(beta, c, mA, alpha_s, cot_beta, Nc), [-1, 0, 1])


class VirtualCorrections:
    """V1 and V2 of the Type-II amplitude.  V1 uses the exact one-loop remainder; V2 uses the
    exact one-loop remainder and the two-loop threshold series (valid for beta_A <~ 0.3)."""

    def __init__(self):
        self.one_loop = OneLoop()
        self.series = ThresholdSeries()

    def at(self, beta, costh, mA=125, muR=125, Nc=3, nf=5):
        """(V1, V2, sum_i |M_i^(0)|^2) at one phase-space point."""
        m0 = [born(i, beta, costh, mA) for i in (1, 2)]
        m1 = [self.one_loop(i, beta, costh, mA, muR, Nc, nf) for i in (1, 2)]
        m2 = [self.series(i, beta, costh, mA, muR, Nc, nf, order=6) for i in (1, 2)]
        norm = sum(abs(a) ** 2 for a in m0)
        v1 = 2 * sum(mp.re(mp.conj(a) * b) for a, b in zip(m0, m1)) / norm
        v2 = (2 * sum(mp.re(mp.conj(a) * b) for a, b in zip(m0, m2)) + sum(abs(b) ** 2 for b in m1)) / norm
        return v1, v2, norm

    def angle_averaged(self, beta, mA=125, muR=125, Nc=3, nf=5, nodes=16):
        """V1 and V2 weighted with the LO angular distribution (Gauss-Legendre in cos(theta))."""
        xs, ws = _gauss_legendre(nodes)
        total = [mp.mpf(0)] * 3
        for c, w in zip(xs, ws):
            v1, v2, norm = self.at(beta, c, mA, muR, Nc, nf)
            total[0] += w * v1 * norm
            total[1] += w * v2 * norm
            total[2] += w * norm
        return total[0] / total[2], total[1] / total[2]


def _gauss_legendre(n):
    """Nodes and weights on [-1, 1]."""
    nodes, weights = [], []
    for k in range(1, n + 1):
        x = mp.cos(mp.pi * (k - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for j in range(2, n + 1):
                p0, p1 = p1, ((2 * j - 1) * x * p1 - (j - 1) * p0) / j
            dp = n * (x * p1 - p0) / (x**2 - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < mp.eps * 10:
                break
        nodes.append(x)
        weights.append(2 / ((1 - x**2) * dp**2))
    return nodes, weights
