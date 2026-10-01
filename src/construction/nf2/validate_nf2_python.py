"""High-precision numerical checks for the Nf^2 coefficient."""

import json
import pathlib
import pickle

import mpmath as mp
import sympy as s


ROOT = pathlib.Path(__file__).parent
x, z, m, mu, nc = s.symbols("x z M2 muR2 nc")
beta, costh, L = s.symbols("beta costh L")

mp.mp.dps = 90

points = [
    ("1/10", "0", 15625, 15625),
    ("1/20", "1/3", 15625, 62500),
    ("1/100", "-2/3", 22500, 5625),
    ("1/5", "2/3", 22500, 22500),
]

rows = []
for name in ("P1", "P2"):
    raw = pickle.loads((ROOT / f"{name}_nf2_expr.pkl").read_bytes())
    compact = pickle.loads((ROOT / f"{name}_nf2_factored.pkl").read_bytes())
    series = pickle.loads((ROOT / f"{name}_nf2_series6.pkl").read_bytes())

    raw_value = s.lambdify((x, z, m, mu, nc), raw, "mpmath", cse=True)
    compact_value = s.lambdify((x, z, m, mu), compact, "mpmath", cse=True)
    series_value = s.lambdify((beta, costh, m, L), series, "mpmath", cse=True)

    for beta_text, costh_text, mass2, mu2 in points:
        beta_value = mp.mpf(beta_text.split("/")[0]) / int(beta_text.split("/")[1])
        if "/" in costh_text:
            num, den = costh_text.split("/")
            costh_value = mp.mpf(num) / int(den)
        else:
            costh_value = mp.mpf(costh_text)

        x_value = (1 - beta_value) / (1 + beta_value)
        z_value = (
            1 + beta_value * beta_value + 2 * beta_value * costh_value
        ) / (1 - beta_value * beta_value)

        full = raw_value(
            x_value,
            z_value,
            mp.mpf(mass2),
            mp.mpf(mu2),
            mp.mpf(3),
        )
        reduced = compact_value(x_value, z_value, mp.mpf(mass2), mp.mpf(mu2))
        threshold = series_value(
            beta_value,
            costh_value,
            mp.mpf(mass2),
            mp.log(mp.mpf(mass2)) - mp.log(mp.mpf(mu2)),
        )

        relative_difference = abs(full - reduced) / max(1, abs(full))
        assert relative_difference < mp.mpf("1e-60"), (
            name,
            beta_text,
            relative_difference,
        )

        row = {
            "projector": name,
            "beta": beta_text,
            "costh": costh_text,
            "mass2": mass2,
            "mu2": mu2,
            "real": mp.nstr(full.real, 75),
            "imag": mp.nstr(full.imag, 75),
            "compact_scaled_difference": mp.nstr(relative_difference, 8),
            "series6_relative_error": mp.nstr(abs(full - threshold) / abs(full), 12),
        }
        rows.append(row)
        print(row, flush=True)

(ROOT / "python_nf2_checks.json").write_text(json.dumps(rows, indent=2))
