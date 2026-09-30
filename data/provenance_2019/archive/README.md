# Archive of the 2019 runs (original calculation of arXiv:1909.08993)

Files from the author's 2019 working directories `Catani ggHH` and `Catani ggAA` (the latter from
`ggjjAA_One-loop_2026.zip`) (the folder name comes from the
gg -> HH template in `template/`, which reads `gg2HHR`; all other files are for gg -> AA).
Large outputs are stored xz-compressed; `SOURCE_SHA256` is the checksum of the original archive.

| path | content | date |
|---|---|---|
| `one_loop/Renorm1.m`, `Renorm2.m` | Z_GG^2 x (tree + mu^-eps x bare one loop), divided by the tree amplitude | Jan 2019 |
| `one_loop/Catani1L1.m`, `Catani1L2.m` | subtraction of 2 I^(1) M^(0) | Jan 2019 |
| `one_loop/Finite1L1.xz`, `Finite1L2.xz` | one-loop remainders as ratios to the tree amplitude | May 2019 |
| `two_loop/Renorm1.m`, `Renorm2.m` | two-loop renormalisation without the factor 64 (P1 includes the O_G-O_J one-loop amplitude with coefficient CJ in the one-loop slot) | Feb-May 2019 |
| `two_loop/Catani2L1.m`, `Catani2L2.m` | two-loop subtraction at x = 1/11, z = 3/11 | Feb 2022 |
| `two_loop/Finite2L1.xz`, `Finite2L2.xz` | two-loop remainders (ratios to the tree amplitude) at that point | Feb 2022 |
| `two_loop/ggjAA11`, `ggjAA12` | one-loop O_G-O_J amplitudes (bare, through O(eps^2)) | May 2019 |
| `two_loop/Renorm1_N.m` | published two-loop renormalisation: 64 Z_GG^2 x amplitude + (8/6) Z_GG Z_GJ x O_G-O_J amplitude; writes `gg2AAR21new` | 19 Jul 2019 |
| `two_loop/Catani2L1_N.m` | published two-loop subtraction; writes the eps^0 and 1/eps coefficients (`diff2m0`, `diff2m1`) | 19 Jul 2019 |
| `two_loop/1byeps_Coef_Num` | record of the 1/eps coefficient with the O_G-O_J counterterm scaled by FLAG, both projections, seven kinematic points | 13 Jun 2019 |

`scripts/check_archive_2019.py` compares these outputs with the published files
(`results/archive_2019.json`): the one-loop outputs carry the mu_R dependence 2 beta0 ln mu_R^2 M^(0)
missing from the published one-loop files and agree with the corrected remainders of
`ggaa.oneloop`; the published two-loop P1 remainder equals 64 x (archive x tree) plus the
O_G-O_J counterterm of the published pipeline, up to 32 C_F times the one-loop O_G-O_J amplitude
(a real multiple of the tree amplitude). This difference is not an uncertainty of the published
remainder: the earlier run's renormalisation constant exceeds Z_GG^2 = Z_as^2 by 2 beta1 a_s^2/eps
and it has no O_G-O_J counterterm, while its infrared subtraction is the same, so its single poles
cannot cancel. In the published version the counterterm was normalised by the cancellation of the
C_F n_f single pole: in `1byeps_Coef_Num` that pole is proportional to (-1 + 96 FLAG) in both
projections at all seven points (the terms proportional to N^3, which do not involve FLAG, are
2.5e-13 and 1.6e-10 of the n_f term at x = 1/11 and grow towards x = 1e-5 in this machine-number
evaluation; they do not affect the value of FLAG), and 64 x 2 x (1/96) = 8/6 is the coefficient in `Renorm1_N.m`
(2 is the coefficient in `../Renorm2L_P1_old.m`). `scripts/check_renorm_2019.py` (`make renorm-2019`)
checks all of this.
