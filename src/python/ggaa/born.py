"""Tree-level Type-II amplitude for g(p1) g(p2) -> A(p3) A(p4).

With the vertex from -(1/8) C_G Phi_A eps_{mu nu rho sigma} G^{a mu nu} G^{a rho sigma},

    V^{ab}_{alpha beta}(k1, k2) = i C_G eps_{mu alpha rho beta} k1^mu k2^rho delta^{ab},

the two diagrams (t- and u-channel gluon exchange) give, after projection onto the
tensors T_1 and T_2 of eqs. (2.7)-(2.8) of the paper and with C_G^2 and delta^{ab} stripped,

    M_1^(0) = i s / 2,
    M_2^(0) = i m_A^2 (x - z)(x^2 + 1)(x z - 1) / (4 x z (x^2 - x z + 1)),

in the variables of arXiv:1909.08993 (s = m_A^2 (1+x)^2/x, u = -m_A^2 z).  This is the
normalization of the loop remainders in the ancillary files: their expansion reads
C_G^2 [M^(0) + a_s M^(1) + a_s^2 M^(2)] delta^{ab}.  Up to phases, the helicity
amplitudes are A(++) = A(--) = C_G^2 M_1 and A(+-) = A(-+) = C_G^2 M_2.

`helicity_amplitudes` evaluates the two diagrams directly with explicit polarization
vectors and is used in the tests as an independent check of the projected formulae.
"""

import itertools

import mpmath as mp
import numpy as np


def kinematics(beta, costh, mA=125):
    b, c = mp.mpmathify(beta), mp.mpmathify(costh)
    x = (1 - b) / (1 + b)
    z = (1 + b**2 + 2 * b * c) / (1 - b**2)
    s = mp.mpf(mA) ** 2 * (1 + x) ** 2 / x
    return x, z, s


def born(i, beta, costh, mA=125):
    """M_i^(0), couplings and colour stripped."""
    x, z, s = kinematics(beta, costh, mA)
    if i == 1:
        return mp.mpc(0, 1) * s / 2
    m2 = mp.mpf(mA) ** 2
    return mp.mpc(0, 1) * m2 * (x - z) * (x**2 + 1) * (x * z - 1) / (4 * x * z * (x**2 - x * z + 1))


_EPS = np.zeros((4, 4, 4, 4))
for _perm in itertools.permutations(range(4)):
    _EPS[_perm] = np.linalg.det(np.eye(4)[list(_perm)])
_METRIC = np.diag([1.0, -1.0, -1.0, -1.0])


def helicity_amplitudes(beta, costh, mA=125.0):
    """Tree amplitudes {(l1, l2): value} from the two Feynman diagrams (double precision,
    C_G^2 and delta^{ab} stripped, overall phase dropped)."""
    beta, costh = float(beta), float(costh)
    s = 4 * mA**2 / (1 - beta**2)
    E, sin = np.sqrt(s) / 2, np.sqrt(1 - costh**2)
    p1, p2 = np.array([E, 0, 0, E]), np.array([E, 0, 0, -E])
    p3 = np.array([E, beta * E * sin, 0, beta * E * costh])
    p4 = p1 + p2 - p3
    low = lambda v: _METRIC @ v

    def exchange(pa):
        q = p1 - pa
        left = np.einsum("kmla,k,l->ma", _EPS, low(p1), low(-q))
        right = np.einsum("knlb,k,l->nb", _EPS, low(p2), low(q))
        return left @ np.linalg.inv(_METRIC) @ right.T / (q @ _METRIC @ q)

    X = exchange(p3) + exchange(p4)
    e1 = {+1: np.array([0, 1, 1j, 0]) / np.sqrt(2), -1: np.array([0, 1, -1j, 0]) / np.sqrt(2)}
    e2 = {+1: np.array([0, 1, -1j, 0]) / np.sqrt(2), -1: np.array([0, 1, 1j, 0]) / np.sqrt(2)}
    return {(l1, l2): e1[l1] @ X @ e2[l2] for l1 in (1, -1) for l2 in (1, -1)}
