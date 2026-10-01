"""Collect the Nf^2 coefficient by powers of the logarithms."""

import pathlib
import pickle
import sys
import time
from functools import lru_cache

import sympy as s


ROOT = pathlib.Path(__file__).parent
name = sys.argv[1]
started = time.time()
expression = pickle.loads((ROOT / f"{name}_nf2_expr.pkl").read_bytes())

mass2 = s.Symbol("M2")
logs = sorted(expression.atoms(s.log), key=str)
dummy_logs = s.symbols("L0:" + str(len(logs)))
log_to_dummy = dict(zip(logs, dummy_logs))
dummy_to_log = dict(zip(dummy_logs, logs))
expression = expression.xreplace(log_to_dummy)

zero_power = (0,) * len(dummy_logs)
log_index = {log: i for i, log in enumerate(dummy_logs)}


@lru_cache(None)
def collect_logs(term):
    if term in log_index:
        powers = list(zero_power)
        powers[log_index[term]] = 1
        return {tuple(powers): s.S.One}

    if not term.has(*dummy_logs):
        return {zero_power: term}

    if term.func == s.Add:
        collected = {}
        for piece in term.args:
            for powers, coefficient in collect_logs(piece).items():
                collected[powers] = collected.get(powers, 0) + coefficient
        return collected

    if term.func == s.Mul:
        collected = {zero_power: s.S.One}
        for piece in term.args:
            next_collected = {}
            for left_power, left_coefficient in collected.items():
                for right_power, right_coefficient in collect_logs(piece).items():
                    powers = tuple(
                        a + b for a, b in zip(left_power, right_power)
                    )
                    next_collected[powers] = (
                        next_collected.get(powers, 0)
                        + left_coefficient * right_coefficient
                    )
            collected = next_collected
        return collected

    if term.func == s.Pow and term.exp.is_Integer and 0 <= term.exp <= 4:
        collected = {zero_power: s.S.One}
        for _ in range(int(term.exp)):
            next_collected = {}
            for left_power, left_coefficient in collected.items():
                for right_power, right_coefficient in collect_logs(term.base).items():
                    powers = tuple(
                        a + b for a, b in zip(left_power, right_power)
                    )
                    next_collected[powers] = (
                        next_collected.get(powers, 0)
                        + left_coefficient * right_coefficient
                    )
            collected = next_collected
        return collected

    raise RuntimeError(str(term)[:300])


log_coefficients = collect_logs(expression)
print(
    name,
    "collected",
    len(log_coefficients),
    time.time() - started,
    flush=True,
)
(ROOT / f"{name}_logcoefs.pkl").write_bytes(
    pickle.dumps((logs, dummy_logs, log_coefficients))
)


@lru_cache(None)
def homogeneous_degree(expression, variable):
    if expression == variable:
        return s.S.One
    if not expression.has(variable):
        return s.S.Zero
    if expression.func == s.Mul:
        return sum(homogeneous_degree(arg, variable) for arg in expression.args)
    if expression.func == s.Pow and not expression.exp.has(variable):
        return expression.exp * homogeneous_degree(expression.base, variable)
    if expression.func == s.Add:
        degrees = {
            homogeneous_degree(arg, variable)
            for arg in expression.args
            if arg != 0
        }
        if len(degrees) == 1:
            return degrees.pop()
    raise ValueError(("nonhomogeneous expression", variable, str(expression)[:160]))


@lru_cache(None)
def normalize_rational(expression):
    if expression.is_Atom:
        return expression
    if expression.func == s.Add:
        return s.cancel(s.Add(*(normalize_rational(arg) for arg in expression.args)))
    if expression.func == s.Mul:
        return s.cancel(s.Mul(*(normalize_rational(arg) for arg in expression.args)))
    if expression.func == s.Pow:
        return s.cancel(normalize_rational(expression.base) ** expression.exp)
    raise ValueError(("unexpected rational function", expression))


terms = []
colour = s.Symbol("nc")
for powers, coefficient in sorted(log_coefficients.items()):
    print(
        "start",
        powers,
        s.count_ops(coefficient),
        time.time() - started,
        flush=True,
    )

    assert homogeneous_degree(coefficient, mass2) == 1
    assert homogeneous_degree(coefficient, colour) == 0

    coefficient = mass2 * normalize_rational(
        coefficient.subs({mass2: 1, colour: 1})
    )
    print(
        "done",
        powers,
        s.count_ops(coefficient),
        coefficient,
        time.time() - started,
        flush=True,
    )

    terms.append(
        coefficient
        * s.prod(log**power for log, power in zip(dummy_logs, powers))
    )

compact = s.Add(*terms).xreplace(dummy_to_log)
(ROOT / f"{name}_nf2_compact.pkl").write_bytes(pickle.dumps(compact))
(ROOT / f"{name}_nf2_compact.txt").write_text(str(compact))
print("FINAL", compact, flush=True)
