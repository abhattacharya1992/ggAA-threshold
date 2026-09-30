# Normalization of the 2019 two-loop remainders

Mathematica steps of the 2019 calculation (projector P1), kept to document the
normalization of the published two-loop files.

| file | step |
|---|---|
| `Renorm2L_P1_old.m` | UV renormalization through two loops; physical C_J; writes `gg2AAR21` |
| `Renorm2L_P1_new.m` | UV renormalization through two loops; multiplies the whole O_G–O_G amplitude by 64 and the O_G–O_J counterterm term by 8/6; writes `gg2AAR21new` |
| `Catani2L_P1_old.m`, `Catani2L_P1_new.m` | infrared subtraction, M^(2),fin = M^(2) − 4 I^(2) M^(0) − 2 I^(1) M^(1), with M^(0) and M^(1) taken from the same renormalized file |

In the `_new` version the tree and one-loop amplitudes used in the infrared subtraction carry
the same factor 64, so the resulting two-loop remainder is 64 times the remainder in the
normalization of the tree amplitudes `data/born/original/ggAA01`, `ggAA02` and of the
published one-loop remainders, which do not have this factor (one-loop step: plain Z_GG²).

The tree normalization is fixed independently: a direct evaluation of the two Feynman
diagrams from −(1/8) C_G Φ_A O_G gives the amplitude C_G² M_i^(0) δ^{ab} with M_i^(0) as in
`ggAA01`, `ggAA02` (`scripts/check_normalization.py`, check 2). The published two-loop
remainders therefore correspond to stripping (C_G/8)², the one-loop remainders and the tree to
stripping C_G². Dividing the two-loop remainders by 64 brings all orders to the same
normalization.

A second, independent defect concerns the published one-loop remainders: their eight
Log[muR2] terms cancel, so they do not depend on μ_R, and under a common rescaling of m_A and
μ_R they change by −2β₀ ln λ² M^(0) (they contain ln m_A² where ln(m_A²/μ_R²) is required).
With the renormalization of `Renorm2L_P1_new.m` (Z_GG² with Z_GG = Z_as, C_G ∝ a_s) and the
I operators of `Catani2L_P1_new.m`, renormalization-group invariance requires

    dM^(1),fin/d ln μ² = 2 β₀ M^(0),
    dM^(2),fin/d ln μ² = 3 β₀ M^(1),fin + (2β₁ − 8 H2g − 12 C_F n_f) M^(0),

where H2g is the constant of `Catani2L_P1_new.m` (no μ factor on that term) and −12 C_F n_f comes
from the (8/6) Z_GG Z_GJ counterterm. The two-loop files divided by 64 satisfy the second
relation, including the M^(0) term, only with

    M^(1),fin(μ_R) = [published] + 2 β₀ ln(μ_R²/GeV²) M^(0),

i.e. the published one-loop files are the remainders at μ_R = 1 GeV
(`scripts/check_normalization.py`, checks 4 and 5; `ggaa.oneloop` applies the term).

The threshold expansion and its validation are unaffected: they compare the series with the
unexpanded two-loop remainder in one and the same normalization.

## Archive of the 2019 runs

`archive/` holds the original Mathematica runs (see `archive/README.md`). `scripts/check_archive_2019.py`
(`make archive`) confirms both defects from them: the archive one-loop remainders carry the term
2 beta0 ln mu_R^2 M^(0) and agree with the corrected remainders to 5e-56, and the published two-loop
P1 remainder equals 64 x (archive x tree) plus the O_G-O_J counterterm, with real parts agreeing to
5e-31. The imaginary parts differ by 32 C_F times the one-loop O_G-O_J amplitude (24 C_F n_f M^(0)
after division by 64): the two versions of the code treat the O_J counterterm differently, which
changes a_s^2 V2 by 0.03 and is part of the open O_J item.
