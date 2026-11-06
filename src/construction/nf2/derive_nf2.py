"""Reduce the Nf^2 terms and construct their beta expansion."""

import json
import pathlib
import pickle

import sympy as s
from sympy.printing.mathematica import mathematica_code


ROOT = pathlib.Path(__file__).parent

x, z, m, mu = s.symbols("x z M2 muR2")
y = (1 + x * x) / x - z
A = s.log(m) - s.log(mu)
B = 2 * s.log(1 + x) - s.log(x)


def H(a, b):
    return (
        216 * a * a
        + 216 * a * b
        + 72 * b * b
        - (480 + 216 * s.I * s.pi) * a
        - (240 + 144 * s.I * s.pi) * b
        + 200
        - 69 * s.pi**2
        + 240 * s.I * s.pi
    )


def J(a, b, t):
    return 36 * t * t + (108 * a + 36 * b - 120 - 36 * s.I * s.pi) * t


r1 = m * (1 + x) ** 2 / x
ry = m * (y * z - 1) / y
rz = m * (y * z - 1) / z

form1 = (
    16
    * s.I
    * r1
    / s.Integer(81)
    * (H(A, B) + J(A, B, s.log(y)) + J(A, B, s.log(z)))
)
form2 = (
    16
    * s.I
    / s.Integer(81)
    * ((ry + rz) / 2 * H(A, B) + ry * J(A, B, s.log(y)) + rz * J(A, B, s.log(z)))
)

checks = []
for name, compact in (("P1", form1), ("P2", form2)):
    full = pickle.loads((ROOT / f"{name}_nf2_compact.pkl").read_bytes())
    logs = sorted((full - compact).atoms(s.log), key=str)
    dummy_logs = s.symbols("L0:" + str(len(logs)))
    difference = s.expand((full - compact).xreplace(dict(zip(logs, dummy_logs))))
    polynomial = s.Poly(difference, *dummy_logs, domain="EX")
    coefficients = [s.cancel(value) for _, value in polynomial.terms()]

    assert all(value == 0 for value in coefficients), (name, coefficients)
    checks.append(
        {
            "name": name,
            "log_coefficient_identities": len(coefficients),
            "all_zero": True,
        }
    )

    (ROOT / f"{name}_nf2_factored.wl").write_text(mathematica_code(compact) + "\n")
    (ROOT / f"{name}_nf2_factored.pkl").write_bytes(pickle.dumps(compact))

print("Exact compact formula comparisons", checks, flush=True)

beta, costh, L = s.symbols("beta costh L")
y_beta = (1 + beta * beta - 2 * beta * costh) / (1 - beta * beta)
z_beta = (1 + beta * beta + 2 * beta * costh) / (1 - beta * beta)
order = 6


def trunc(expression, power=order):
    return s.series(expression, beta, 0, power + 1).removeO().expand()


B_beta = s.log(4) + sum(
    beta ** (2 * k) / s.Integer(k) for k in range(1, order // 2 + 1)
)
log_y = trunc(s.log(y_beta))
log_z = trunc(s.log(z_beta))
H_beta = trunc(H(L, B_beta))
J_y = trunc(J(L, B_beta, log_y))
J_z = trunc(J(L, B_beta, log_z))
r_y = trunc(m * (y_beta * z_beta - 1) / y_beta)
r_z = trunc(m * (y_beta * z_beta - 1) / z_beta)

series1 = trunc(
    64
    * s.I
    * m
    / s.Integer(81)
    * sum(beta ** (2 * k) for k in range(order // 2 + 1))
    * (H_beta + J_y + J_z)
)
series2 = trunc(
    16
    * s.I
    / s.Integer(81)
    * ((r_y + r_z) * H_beta / 2 + r_y * J_y + r_z * J_z)
)

for name, series in (("P1", series1), ("P2", series2)):
    coefficients = {
        power: s.factor(series.coeff(beta, power)) for power in range(order + 1)
    }
    print(name, coefficients, flush=True)

    output = sum(value * beta**power for power, value in coefficients.items())
    (ROOT / f"{name}_nf2_series6.wl").write_text(mathematica_code(output) + "\n")
    (ROOT / f"{name}_nf2_series6.pkl").write_bytes(pickle.dumps(output))
    (ROOT / f"{name}_nf2_series6.tex").write_text(s.latex(output))

(ROOT / "nf2_exact_reduction_checks.json").write_text(json.dumps(checks, indent=2))
