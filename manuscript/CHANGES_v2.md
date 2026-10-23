# Changes from arXiv v1 to v2 (arXiv:2610.19024)

The arXiv v1 source is the tag `arxiv-v1` (commit 097cce7) of this repository; its
`manuscript/main.pdf` is the arXiv v1 PDF. Compare with `git diff arxiv-v1 -- manuscript`.

The notation changes make each symbol stand for one quantity only.
No change to the numerical data is intended by the notation changes below.

1. **The variable of the amplitude files is now v, not x.** Eqs. (2.2), (2.4), (2.12), Sec. 3.1
   (the map and the Li₂ example) and App. B (`v = 1/11`); a half-sentence after eq. (2.2) says that
   v is the x of arXiv:1909.08993, renamed because x is the momentum fraction of Sec. 4. The
   sentence "Here x and z are momentum fractions, not the variables of Eq. (2.2)" is gone.
2. **The argument of the g → gg factors is ξ** (eq. 4.8 and the two sentences around it), so that
   z is only the file variable.
3. **The n_f² threshold constant of eqs. (3.6)–(3.7) is K_{n_f}** (it was K, which also denotes
   the one-loop forward limit of M₂).
4. **Eq. (4.17) names the function: W(x) = [x⁴ + (1−x)⁴]/[x(1−x)]**, and the sentence after it
   gives W(x) in terms of the factors of eq. (4.8).
5. **"Weights" are no longer called weights.** The w's of eq. (4.8) are "the longitudinal factors
   of the light cone wavefunctions" ("wavefunction factors" for short; paragraph title "Helicity
   selection and wavefunction factors"); F_A and F_S are "helicity sums", defined after eq. (4.10)
   as sums over daughter helicities of products of two wavefunction factors, overlaps of amplitudes
   rather than probabilities. Changed in the abstract, Sec. 1 (item 2), Sec. 4.3 (four places), the
   caption of Fig. 4, Sec. 4.4, Sec. 6 and App. D. The label `eq:splitting-weights` is unchanged.
6. **ℒ and 𝒢.** After eq. (4.2): "Away from threshold the natural variable is
   𝒢 = ½ ln[(m_A⁴ − tu)/m_A⁴] = ℒ − ln(1 − β_A²), used in Sec. 4.4." (Sec. 4.4 already had the
   relation.)

The package code is untouched: `src/python/ggaa` keeps `x` as the name of the file variable
internally, as the published files do.

## Follow-up corrections on 22 October 2026

- Corrected a remaining amplitude-file variable from `x` to `v` and used the
  defined term "helicity sums" consistently in the current physics section.
- Described the nonresonant ALP result cited from Brivio et al. as an expected
  sensitivity, rather than an observed bound, and added the Borowka et al.
  erratum to the references.
- Restored the qualifiers of the arXiv v1 abstract.
- Labelled the M₁ leading-power line, eq. (4.3), and cited (4.3)–(4.4).
- Added Manohar–Waalewijn, Gaunt (2013) and Buffing–Diehl–Kasemets (with erratum)
  to the SPS/DPS overlap sentence.
