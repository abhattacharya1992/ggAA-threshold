# Archive of the 2019 runs (original calculation of arXiv:1909.08993)

Files from the author's 2019 working directory `Catani ggHH` (the folder name comes from the
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

`scripts/check_archive_2019.py` compares these outputs with the published files
(`results/archive_2019.json`): the one-loop outputs carry the mu_R dependence 2 beta0 ln mu_R^2 M^(0)
missing from the published one-loop files and agree with the corrected remainders of
`ggaa.oneloop`; the published two-loop P1 remainder equals 64 x (archive x tree) plus the
O_G-O_J counterterm of the published pipeline, up to 32 C_F times the one-loop O_G-O_J amplitude
(a real multiple of the tree amplitude).
