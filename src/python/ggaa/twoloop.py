"""Numerical evaluation of the published two-loop Type-II finite remainders.

The ancillary files Finite_2Loop_P1 and Finite_2Loop_P2 of arXiv:1909.08993 are single
Wolfram Language expressions of 70 and 133 MB in x, z, Y = 1 + y, Z = 1 + z,
M2 = m_A^2, muR2 = mu_R^2, N = N_c and Nf, built from Log, PolyLog, Zeta and Li22.
They are evaluated here in one pass over the text: an operator-precedence (shunting-yard)
evaluator reads the tokens as they come and keeps only a stack of numbers, so memory
use stays small.  Special functions are evaluated with mpmath, and Li22 with
`ggaa.li22`, and repeated arguments are cached.

This code shares nothing with the Wolfram Language implementation used for the
validation table; agreement between the two is an independent check of the reference
values.
"""

import json
import re

import mpmath as mp

from .li22 import li22
from .paths import INPUTS

FILES = {(1, 1): "Finite_1Loop_P1", (2, 1): "Finite_1Loop_P2",
         (1, 2): "Finite_2Loop_P1", (2, 2): "Finite_2Loop_P2"}

_TOKEN = re.compile(r"\s*(?:(\d+)|([A-Za-z][A-Za-z0-9]*)(\[)?|(.))")
_PRECEDENCE = {"+": 1, "-": 1, "*": 2, "/": 2, "neg": 3, "^": 4}


def _apply(op, values):
    if op == "neg":
        values.append(-values.pop())
        return
    b, a = values.pop(), values.pop()
    if op == "+":
        values.append(a + b)
    elif op == "-":
        values.append(a - b)
    elif op == "*":
        values.append(a * b)
    elif op == "/":
        values.append(a / b)
    else:
        values.append(a ** b if not (isinstance(b, mp.mpf) and b == int(b)) else a ** int(b))


class Evaluator:
    """Evaluate one Wolfram Language expression for given symbol values.

    If `store` is a path, computed Li22 values are also kept there, so that a long
    evaluation interrupted by the user can resume without recomputing them."""

    def __init__(self, symbols, store=None):
        self.symbols = symbols
        self._cache = {}
        self._store_path = store
        self._store = json.loads(store.read_text()) if store and store.exists() else {}

    def function(self, name, args):
        key = (name,) + tuple(mp.nstr(a, mp.mp.dps - 8) for a in args)
        if key not in self._cache:
            if name == "Log":
                value = mp.log(args[0])
            elif name == "PolyLog":
                value = mp.polylog(int(args[0]), args[1])
            elif name == "Zeta":
                value = mp.zeta(int(args[0]))
            elif name == "Li22":
                label = f"{mp.mp.dps}:" + ",".join(mp.nstr(a, mp.mp.dps - 8) for a in args)
                if label in self._store:
                    value = mp.mpmathify(self._store[label])
                else:
                    value = li22(*args)
                    if self._store_path:
                        self._store[label] = mp.nstr(value, mp.mp.dps)
                        self._store_path.write_text(json.dumps(self._store))
            else:
                raise ValueError(f"unknown function {name}")
            self._cache[key] = value
        return self._cache[key]

    def evaluate(self, text):
        values, ops, arity = [], [], []
        expect_operand = True
        for match in _TOKEN.finditer(text):
            number, name, bracket, other = match.groups()
            if number:
                values.append(mp.mpf(number))
                expect_operand = False
            elif name and bracket:
                ops.append(("call", name))
                arity.append(1)
                expect_operand = True
            elif name:
                values.append(self.symbols[name])
                expect_operand = False
            elif other in "([":
                ops.append(other)
                expect_operand = True
            elif other == ")":
                while ops[-1] != "(":
                    _apply(ops.pop(), values)
                ops.pop()
                expect_operand = False
            elif other == ",":
                while not isinstance(ops[-1], tuple):
                    _apply(ops.pop(), values)
                arity[-1] += 1
                expect_operand = True
            elif other == "]":
                while not isinstance(ops[-1], tuple):
                    _apply(ops.pop(), values)
                _, fname = ops.pop()
                n = arity.pop()
                args = values[-n:]
                del values[-n:]
                values.append(self.function(fname, args))
                expect_operand = False
            elif other in "+-*/^":
                if expect_operand:
                    if other == "-":
                        ops.append("neg")
                    continue
                prec = _PRECEDENCE[other]
                while ops and isinstance(ops[-1], str) and ops[-1] in _PRECEDENCE:
                    top = _PRECEDENCE[ops[-1]]
                    if top > prec or (top == prec and other != "^"):
                        _apply(ops.pop(), values)
                    else:
                        break
                ops.append(other)
                expect_operand = True
            elif other:
                raise ValueError(f"unexpected character {other!r}")
        while ops:
            _apply(ops.pop(), values)
        if len(values) != 1:
            raise ValueError("malformed expression")
        return values[0]


def read_expression(path):
    """Right-hand side of a published file; long integers split with a trailing backslash are rejoined."""
    return path.read_text().split("=", 1)[1].replace("\\\n", "")


def symbols(beta, costh, mA=125, muR=125, Nc=3, nf=5):
    b, c = mp.mpmathify(beta), mp.mpmathify(costh)
    x = (1 - b) / (1 + b)
    z = (1 + b**2 + 2 * b * c) / (1 - b**2)
    y = x + 1 / x - z
    return {"x": x, "z": z, "Y": 1 + y, "Z": 1 + z, "M2": mp.mpf(mA) ** 2, "muR2": mp.mpf(muR) ** 2,
            "N": mp.mpf(Nc), "Nf": mp.mpf(nf), "Pi": +mp.pi, "I": mp.mpc(0, 1)}


def remainder(i, loops, beta, costh, mA=125, muR=125, Nc=3, nf=5, directory=INPUTS, store=None):
    """M_i^{II,(loops),fin} from the published file, couplings stripped, in the normalization
    of that file.  For loops = 2 this is 64 times the package normalization (ggaa.normalization)."""
    text = read_expression(directory / FILES[i, loops])
    return Evaluator(symbols(beta, costh, mA, muR, Nc, nf), store).evaluate(text)
