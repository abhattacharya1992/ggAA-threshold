# Outputs of the construction and validation runs

`data/construction` and `data/validation` were produced in September 2026 by Mathematica
runs with the code in `src/construction`. This folder keeps the outputs of those runs,
compressed with xz.

| file | run | content |
|---|---|---|
| `known/P1_nf2_p4.wl`, `known/P1_nf2_p6.wl`, `known/P2_nf2_p4.wl`, `known/P2_nf2_p6.wl` | `nf2/` | coefficients of n_f² in M1 and M2 through β⁴ and β⁶ |
| `known/P1_nf1_p4.wl`, `known/P1_nf1_p6.wl` | `p1_nf1/` | coefficient of n_f in M1 through β⁴ and β⁶ |
| `P1_nf0_series4.wl`, `P1_nf0_series6.wl`, `P2_nf0_series4.wl`, `P2_nf0_series6.wl`, `P2_nf1_series4.wl`, `P2_nf1_series6.wl` | `remaining/` | coefficients of n_f⁰ in M1 and M2 and of n_f in M2 |
| `validation_all.csv` | `validation/` | the 1080 validation records |

All files are in the normalisation of the published remainders. `SHA256SUMS` gives the SHA-256
of the uncompressed files.

`make provenance` checks the outputs against `SHA256SUMS`, checks that every file in
`data/construction` agrees exactly with the corresponding coefficient of these outputs after
both are brought to the same form, and checks that the 1080 records in
`data/validation/validation_1080.csv` agree with `validation_all.csv` to the 40 digits stored.
