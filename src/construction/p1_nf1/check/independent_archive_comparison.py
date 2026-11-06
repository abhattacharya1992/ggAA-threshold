"""Exact comparison of the saved Wolfram series with the Python recurrence."""

import argparse
import json
import re
from pathlib import Path

import sympy as s
from sympy.parsing.mathematica import parse_mathematica


def split_top(text, separator=","):
    """Split at top-level separators, leaving nested Wolfram syntax intact."""
    stack = []
    pieces = []
    start = 0
    quoted = False
    escaped = False
    i = 0

    while i < len(text):
        char = text[i]

        if quoted:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == '"':
                quoted = False
            i += 1
            continue

        if char == '"':
            quoted = True
            i += 1
            continue

        if text[i : i + 2] == "<|":
            stack.append("|>")
            i += 2
            continue
        if text[i : i + 2] == "|>":
            assert stack.pop() == "|>"
            i += 2
            continue

        if char in "[{(":
            stack.append({"[": "]", "{": "}", "(": ")"}[char])
        elif char in "]})":
            assert stack.pop() == char
        elif text[i : i + len(separator)] == separator and not stack:
            pieces.append(text[start:i].strip())
            start = i + len(separator)
            i += len(separator)
            continue

        i += 1

    assert not stack and not quoted
    pieces.append(text[start:].strip())
    return pieces


def parse_association(text):
    text = text.strip()
    assert text.startswith("<|") and text.endswith("|>")

    result = {}
    for item in split_top(text[2:-2]):
        if not item:
            continue
        key_value = split_top(item, "->")
        assert len(key_value) == 2
        key, value = key_value
        key = key[1:-1] if key.startswith('"') else int(key)
        result[key] = value
    return result


def parse_expression(text):
    # Avoid SymPy interpreting beta while the Mathematica expression is parsed.
    expression = parse_mathematica(
        re.sub(r"\bbeta\b", "betavar", text.replace("\\\n", ""))
    )
    return expression.xreplace({s.Symbol("betavar"): s.Symbol("beta")})


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Exact rational comparison of saved Wolfram coefficients "
            "with their reference."
        )
    )
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    root = args.run_directory.resolve()
    reference_record = parse_association(
        (root / "source/python_reference_series.wl").read_text()
    )
    references = {
        key: parse_expression(value) for key, value in reference_record.items()
    }

    beta = s.Symbol("beta")
    rows = []

    for order in (4, 6):
        record = parse_association(
            (root / f"results/P1_nf1_series{order}.wl").read_text()
        )
        actual = {
            key: parse_expression(value)
            for key, value in parse_association(record["NormalizedCoefficients"]).items()
        }
        expected = s.Poly(references[f"NormalizedOrder{order}"], beta)

        for power in range(order + 1):
            difference = s.cancel(actual.get(power, s.S.Zero) - expected.nth(power))
            rows.append(
                {
                    "order": order,
                    "power": power,
                    "difference": str(difference),
                    "pass": bool(difference == 0),
                }
            )

    report = {
        "checks": rows,
        "all_pass": all(row["pass"] for row in rows),
    }
    text = json.dumps(report, indent=2)
    print(text)

    if args.report:
        args.report.write_text(text + "\n")

    return 0 if report["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
