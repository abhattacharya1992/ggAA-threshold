# One-loop gg → HH through two CP-even effective vertices: the forward logarithm

Test of Sec. 4.3: for a CP-even coupling, H G^a_{μν} G^{aμν}, the double parton scattering configuration
gives the scalar helicity weight F_S(x) = 1/[x(1−x)] + 2x(1−x) instead of F_A, x = (1+β)/2, with the
same constant as for AA once each coefficient is expressed in units of its tree amplitude at threshold.

## How the numbers were obtained

The one-loop Type II amplitudes for gg → HH at m_H = m_A were computed by the author with the one-loop
gg → AA code of arXiv:1909.08993 (QGRAF diagrams, FORM projections, LiteRed reduction, Mathematica),
after replacing the effective vertices by those of −(1/4) C Φ G^a_{μν} G^{aμν} (Φgg and Φggg).
Versions: `VERSIONS.txt`. Conventions: m = μ_R = 1, N_c = 3, n_f = 5, couplings set to 1.

Checks recorded in `hh_dps.log`:
- [4] the same code reproduces the archived 2019 one-loop gg → AA result to 40 digits (this fixes the
  branch of √(−1/ε²) in the master integrals; the other branch fails at O(1));
- [5] the HH double poles equal those of Catani's operator for both projections at two angles (15 digits);
  HH tree at threshold M_1^(0) = −2i m² (AA: +2i m²);
- [6] forward-limit coefficient of ln p_T², from the finite parts at 1 − c = 10⁻¹⁴ and 10⁻¹⁸.
  In block [6] the slopes are normalised by their values at β = 10⁻⁶, where F_A = 3.5 + 4.5β² and
  F_S = 4.5 + 3.5β²; the resulting offsets of about 10⁻¹² are these β² terms, not deviations.

## Files

- `dps_hh.csv`: β; forward coefficient of ln p_T² divided by π N_c m² for AA and HH; F_A, F_S at
  x = (1+β)/2; AA/F_A − 1; −HH/F_S − 1 (the sign is that of the HH tree); |M_2 slope|/|M_1 slope| for HH.
  The deviations of order 10⁻¹⁵ come from the finite-difference step (ln(1+c) changes by ~10⁻¹⁴).
- `hh_dps.log`: the full output of the analysis script.
- `VERSIONS.txt`: program versions.

Not included in the calculation: the diagram with the four-gluon term of the CP-even operator (a bubble
carrying the momentum of one Higgs boson) and the contact counterterm of two such vertices. Neither
depends on p_T, so neither affects the logarithm; both are needed only for the complete finite remainder.

## Code

The FORM and Mathematica code of this calculation, which builds on the one-loop code of
arXiv:1909.08993, is not part of this package. It is available from the author on request.

`scripts/check_dps_hh.py` (`make dps-hh`) reads these files and writes `results/dps_hh.json`, which
supplies the numbers of Sec. 4.3 and the HH points of Fig. 4.
