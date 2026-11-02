#!/usr/bin/env python3
"""Real part of the one-loop remainder M_1 in closed form, Eq. (4.12) of the paper.

The published file Finite_1Loop_P1 of arXiv:1909.08993 (run inputs/fetch_inputs.sh first)
depends on x, z, Y = 1 + y and Z = 1 + z, with y = (1 + x^2)/x - z.  On the threshold map
of ggaa.oneloop,

    z - x = 2b(1+c)/(1-b^2),   y - x = 2b(1-c)/(1-b^2),   1 - xz = 2b(1-c)/(1+b)^2,
    1 - xy = 2b(1+c)/(1+b)^2,

so for 0 < b < 1 and -1 < c < 1 every logarithm in the file has a positive argument and
every dilogarithm, Li2(x/z), Li2(xz), Li2(1/Y), Li2(1/Z), Li2(x/y), Li2(xy), an argument
between 0 and 1.  All these functions are real, and the real part of the remainder follows
from the explicit factors of I.  (The scale term added in ggaa.oneloop is proportional to
M_1^(0) = i s/2 and does not contribute.)

The expression is evaluated as a polynomial in the logarithms of the irreducible factors,
the dilogarithms, pi and I (with I^2 = -1), with coefficients in the field
Q(x, z, M2, muR2, N, Nf).  Arithmetic in this field is exact, so the difference between the
real part and Eq. (4.12) vanishes if and only if it is identically zero.
"""

import ast
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

from sympy import QQ
from sympy.polys.fields import field

from ggaa.paths import INPUTS

K, x, z, M2, muR2, N, Nf = field("x,z,M2,muR2,N,Nf", QQ)
RATIONAL = {"x": x, "z": z, "Y": (1 + x + x**2 - x * z) / x, "Z": 1 + z,
            "M2": M2, "muR2": muR2, "N": N, "Nf": Nf}

# Functions of the file in terms of logarithms of irreducible factors; all arguments are
# positive in the physical region.
FUNCTIONS = [
    ("PolyLog[2, x/z]", "Li2_x_over_z"), ("PolyLog[2, x*z]", "Li2_xz"),
    ("PolyLog[2, Y^(-1)]", "Li2_1_over_Y"), ("PolyLog[2, Z^(-1)]", "Li2_1_over_Z"),
    ("PolyLog[2, x/((1 + x^2)/x - z)]", "Li2_x_over_y"), ("PolyLog[2, x*((1 + x^2)/x - z)]", "Li2_xy"),
    ("Log[-x + (1 + x^2)/x - z]", "(ln_1_minus_xz - ln_x)"),       # y - x = (1 - xz)/x
    ("Log[(1 + x^2)/x - z]", "(ln_1_plus_x2_minus_xz - ln_x)"),    # y = (1 + x^2 - xz)/x
    ("Log[Y]", "(ln_1_plus_x_plus_x2_minus_xz - ln_x)"),           # Y = (1 + x + x^2 - xz)/x
    ("Log[-x + z]", "ln_z_minus_x"), ("Log[1 + x]", "ln_1_plus_x"), ("Log[x]", "ln_x"),
    ("Log[z]", "ln_z"), ("Log[Z]", "ln_1_plus_z"), ("Log[M2]", "ln_M2"), ("Log[muR2]", "ln_muR2"),
]
SYMBOLS = {name for _, text in FUNCTIONS for name in re.findall(r"[A-Za-z_][\w]*", text)} | {"PI", "I"}


def multiply(a, b):
    """Product of two monomials (sorted tuples of (symbol, power)) and its sign from I^2 = -1."""
    powers = dict(a)
    for name, p in b:
        powers[name] = powers.get(name, 0) + p
    sign = 1
    if powers.get("I", 0) >= 2:
        q, powers["I"] = divmod(powers["I"], 2)
        sign = -1 if q % 2 else 1
    return tuple(sorted((k, p) for k, p in powers.items() if p)), sign


class Poly(dict):
    """Polynomial in the transcendental symbols, {monomial: coefficient in K}."""

    def __add__(self, other):
        out = Poly(self)
        for k, c in other.items():
            out[k] = out.get(k, K.zero) + c
        return out

    def __neg__(self):
        return Poly({k: -c for k, c in self.items()})

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        out = Poly()
        for k1, c1 in self.items():
            for k2, c2 in other.items():
                k, sign = multiply(k1, k2)
                out[k] = out.get(k, K.zero) + sign * c1 * c2
        return out

    def scalar(self):
        assert set(self) <= {()}, "transcendental function in a denominator"
        return self.get((), K.zero)


def constant(c):
    return Poly({(): K(c) if isinstance(c, int) else c})


def evaluate(node):
    if isinstance(node, ast.Expression):
        return evaluate(node.body)
    if isinstance(node, ast.Constant):
        assert isinstance(node.value, int), "non-integer constant"
        return constant(node.value)
    if isinstance(node, ast.Name):
        return Poly({((node.id, 1),): K.one}) if node.id in SYMBOLS else constant(RATIONAL[node.id])
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        value = evaluate(node.operand)
        return -value if isinstance(node.op, ast.USub) else value
    if isinstance(node, ast.BinOp):
        left, right = evaluate(node.left), evaluate(node.right)
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Div):
            d = right.scalar()
            return Poly({k: c / d for k, c in left.items()})
        if isinstance(node.op, ast.Pow):
            e = right.scalar()
            n = int(e.numer.LC) if e.denom == 1 and e.numer.is_ground else None
            assert n is not None and K(n) == e, "non-integer power"
            if set(left) <= {()}:
                return constant(left.scalar() ** n)
            assert n > 0, "negative power of a transcendental function"
            out = constant(1)
            for _ in range(n):
                out = out * left
            return out
    raise ValueError(f"unexpected expression: {ast.dump(node)[:200]}")


def published(path=INPUTS / "Finite_1Loop_P1"):
    """The published remainder as a Poly."""
    text = path.read_text().split("=", 1)[1].strip().rstrip(";")
    for old, new in FUNCTIONS:
        text = text.replace(old, new)
    assert "Log[" not in text, "unmapped function in the published file"
    text = text.replace("[", "(").replace("]", ")").replace("^", "**")
    text = re.sub(r"\bPi\b", "PI", text)
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 100000))
    value = evaluate(ast.parse("(" + text + ")", mode="eval"))
    return Poly({k: c for k, c in value.items() if c != 0})


def real_part(value):
    return Poly({k: c for k, c in value.items() if dict(k).get("I", 0) == 0})


def closed_form():
    """Eq. (4.12) in the file variables: 1 - b^2 = 4x/(1+x)^2 and
    b^2(1-c^2)/(1-b^2 c^2) = (z-x)(y-x)/(Y Z)."""
    f_a = (1 + x) ** 2 / x - 2 * x / (1 + x) ** 2                      # F_A((1+b)/2)
    log_ratio = Poly({(("ln_z_minus_x", 1),): K.one, (("ln_1_minus_xz", 1),): K.one,
                      (("ln_1_plus_x_plus_x2_minus_xz", 1),): -K.one, (("ln_1_plus_z", 1),): -K.one})
    beta0_term = M2 * (11 * N - 2 * Nf) / 3 * (1 + x) ** 2 / (2 * x)    # m^2 * 2 beta0/(1-b^2)
    return Poly({(("PI", 1),): K.one}) * (constant(M2 * N * f_a) * log_ratio - constant(beta0_term))


def difference():
    """Nonzero terms of Re M_1^(1) - Eq. (4.12); empty if the identity holds."""
    diff = real_part(published()) - closed_form()
    return {k: c for k, c in diff.items() if c != 0}


def main():
    value = published()
    real = real_part(value)
    print(f"Published one-loop remainder M_1: {len(value)} terms in the logarithms, dilogarithms, pi and I.")
    print("Real part (all functions real in the physical region):")
    for k in sorted(real, key=str):
        name = " * ".join(f"{s}^{p}" if p > 1 else s for s, p in k) or "1"
        print(f"   {name}:  {real[k]}")
    diff = difference()
    print("\nRe M_1^(1) - Eq. (4.12):", "identically zero" if not diff else f"nonzero terms {sorted(diff, key=str)}")
    if diff:
        sys.exit(1)


if __name__ == "__main__":
    main()
