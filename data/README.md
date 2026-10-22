# Data

## construction/

`M{i}_nf{k}_beta_{n}.wl` is the coefficient of β_A^n in the coefficient of n_f^k of
M_i^{II,(2),fin}, as produced by the threshold construction (n = m4 … m1 denote
n = −4 … −1). Files are Wolfram Language expressions in `beta` (or `Lb` = ln β_A),
`costh`, `m2`, `mu2`, `nc`, and contain the endpoint constants `MB[p,q,a,w]`.
Inverse and odd powers are present in this form; `scripts/reduce_series.py` shows that
all of them vanish once the constants are reduced. The files were produced by the code in
`src/construction`; `make provenance` checks them against the outputs of those runs in `runs/`.

## runs/

Outputs of the construction and validation runs, compressed with xz; see `runs/README.md`.

## series/

`M{i}_nf{k}.m` is the final series: a Wolfram Language list `{c0, c2, c4, c6}` with

    M_{i,k} = c0 + beta^2 c2 + beta^4 c4 + beta^6 c6 + O(beta^8),

each c_n a polynomial in `Log[beta]`. The complete coefficient is
M_i = M_{i,0} + n_f M_{i,1} + n_f^2 M_{i,2}. The files are read by
`src/python/ggaa/series.py` and `src/mathematica/evaluate_series.wl`, and are written
by `scripts/reduce_series.py`.

## validation/

`validation_1080.csv`: comparison of the series with the unexpanded remainders.

| column | meaning |
|---|---|
| `quantity` | `M{i}_nf{k}` for the coefficient of n_f^k, `M{i}_nf5` for the complete n_f = 5 sum |
| `mA_GeV`, `muR_GeV`, `Nc` | (125, 125, 3), (125, 250, 3) or (150, 75, 4) |
| `beta`, `costheta` | exact rationals: β_A ∈ {1/1000, 1/500, 1/200, 1/100, 1/50, 1/20, 1/10, 1/5, 3/10}, cos θ ∈ {0, ±1/3, ±2/3} |
| `re_exact`, `im_exact` | unexpanded remainder, evaluated at increasing working precision until stable; an independent evaluation at 160 digits (`make exact-check`) reproduces the stored values to 1.6×10⁻⁴⁰, i.e. to the last stored digit |
| `re_series4`, `im_series4` | series through β_A⁴ |
| `re_series6`, `im_series6` | series through β_A⁶ |

Values are given to 40 significant digits. 135 points × 8 quantities = 1080 rows.

`forward_angles.csv`: the same comparison towards the beam direction, for the complete
n_f = 5 coefficients at m_A = μ_R = 125 GeV, N_c = 3: β_A ∈ {1/10, 1/5, 3/10} and
cos θ ∈ {9/10, 99/100, 999/1000, −99/100}. Relative errors were computed at high precision;
the amplitude values are given in double precision. The rows at cos θ = ±99/100 agree, as
required by Bose symmetry. All errors are below the maximum on the main grid at the same β_A, except
for M1 at β_A = 1/10, where the value at cos θ = 999/1000 exceeds it by 21% (Sec. 3.3).

`one_loop_reference.csv`: the published one-loop remainders at six points, evaluated at
75 digits, used by `scripts/check_one_loop.py` and the tests.

## dps_hh/

One-loop gg → HH through two CP-even effective vertices: forward coefficient of ln p_T² against the scalar
weight F_S (Sec. 4.3, Fig. 4). How the numbers were obtained, the checks and the files: `dps_hh/README.md`.

## constants/

`mb_constants.csv`: the thirteen constants MB[p,q,a,w] of Appendix A to 60 digits and
their reductions to ζ₃, π, ln 2 and Li₄(½). Written by `scripts/check_mb_constants.py`.

**Normalization.** `validation_1080.csv` and `forward_angles.csv` store two-loop values in the
normalization of the published 2019 files, which is 64 times the package normalization
(`src/python/ggaa/normalization.py`); `ggaa.validation.records` divides by 64 on reading.
Relative errors do not depend on this factor. `data/construction` is also in the published
normalization; `data/series` is in the package normalization.
