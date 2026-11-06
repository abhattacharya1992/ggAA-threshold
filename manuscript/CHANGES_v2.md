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

## Corrections on 2 November 2026

Equation numbers refer to the corrected version. Eq. (4.12) is new, so the later equations
of Sec. 4 are renumbered by one.

- **p_T dependence of the omitted HH diagrams** (Sec. 4.3, after eq. (4.13), and
  `data/dps_hh/README.md`). The diagram with the four-gluon term of the CP-even operator
  does depend on p_T: its bubble gives a term proportional to (ε₁·q)(ε₂·q) ∝ p_T², which is
  regular at p_T = 0. "Neither depends on p_T" is replaced by "Neither produces a logarithm
  of p_T as p_T → 0". The tested logarithmic coefficient is not affected.
- **Cross-references.** Sec. 1: V₂ is defined in eq. (5.2) (the text said "Sec. 5.3").
  App. B: I⁽²⁾ and H_g⁽²⁾ are eqs. (3.3) and (3.5) of arXiv:1909.08993 (the text cited
  eq. (3.4), which defines K).
- **Errata in the bibliography.** JHEP.bst does not print the `note` field, so the errata of
  Dixon–Esterlis (JHEP 08 (2016) 131), Diehl–Ostermeier–Schäfer (JHEP 03 (2016) 001) and
  Plehn–Spira–Zerwas (Nucl. Phys. B 531 (1998) 655) did not appear in the bibliography. They
  are now separate entries, cited with the original papers, as for Borowka et al. and
  Buffing et al.
- **Status of the endpoint relations** (Sec. 3.2, Sec. 6). The relations between the endpoint
  constants are identified numerically with PSLQ and have not been proved analytically.
  Sec. 3.2 now says so and describes the cancellations as exact once these relations are
  used; its title is "Endpoint constants and cancellations".
- **Data availability.** "form a self-contained package" is replaced by "form one package, in
  which the one-loop gg → HH coefficients enter as computed values". The code that generated
  these coefficients is available on request, as before.
- **Closed form of Re M₁⁽¹⁾** (Sec. 4.3, new eq. (4.12)):
  Re M₁⁽¹⁾ = π m² [N_c F_A((1+β)/2) ln(β²(1−c²)/(1−β²c²)) − 2β₀/(1−β²)].
  It agrees with the published one-loop remainder at 1350 points with 0.001 ≤ β ≤ 0.99,
  |c| ≤ 0.999, (N_c, n_f) = (3,5), (3,0), (4,2), (5,6), (2,1) and
  (m_A, μ_R) = (125,125), (125,250), (300,40) GeV; with 120-digit arithmetic the relative
  differences are below 10⁻¹⁰⁵. Code: `ggaa.oneloop.real_m1`, part 4 of
  `scripts/check_one_loop.py` (results in `results/one_loop.json`), and the test
  `OneLoop.test_real_part_closed_form`. The closed form implies eq. (4.11) for all β and
  gives Re Â₁₀⁽¹⁾ = −π(7N_c ln 2 + 2β₀).
- **v₁ and ℓ_max** (Sec. 5.2). With this constant the β₀ terms cancel, so
  v₁ = −(49/2)π²N_c² ln 2, independent of n_f and μ_R, and ℓ_max = −7 ln 2, that is,
  p_T = m_A/256. The numerical values are unchanged.
- **Fig. 5(b)** (file `fig6_cross_section.pdf`). The axis label and the legend show the
  angular averages a_s⟨V₁⟩ and a_s²⟨V₂⟩, as in the caption (`scripts/make_physics_figures.py`).

## Further corrections on 2 November 2026

The acknowledgments are not changed here.

- **Abstract.** ln³p_T and ln²p_T are written as ln³(2p_T/m_A) and
  ln²(2p_T/m_A), and the sentence on planar N = 4 super Yang–Mills theory states that the
  ratio of the leading two- and one-loop logarithms agrees. The forward-limit agreement of the
  one-loop coefficient is stated as exact (eq. (4.12)) instead of "to 3×10⁻¹⁵".
- **Eq. (4.7).** The corrections are written as O_log(β_A²) instead of O(β_A). A new sentence
  states that they are of order β_A² ln β_A for M₁⁽¹⁾ and β_A² for M₂⁽¹⁾, from evaluations at
  10⁻⁶ ≤ β_A ≤ 10⁻³ (part 5 of `scripts/check_one_loop.py`).
- **Eq. (4.12) is derived.** In the physical region all logarithms of the published one-loop
  remainder have positive arguments and all dilogarithms have arguments between 0 and 1, so the
  real part follows from the explicit factors of i. `scripts/derive_one_loop_real.py` evaluates
  it in exact arithmetic over Q(x, z, m_A², μ_R², N_c, n_f): it contains no dilogarithms and
  equals eq. (4.12) identically (output in `results/one_loop_real.txt`, run by `make one-loop`;
  test `OneLoop.test_real_part_derivation`). Sec. 4.3 now says so and describes the comparison
  at 1350 points as a numerical check.
- **Eq. (4.18).** "The forward limit of M₂⁽¹⁾ itself has a closed form. The exact one-loop
  remainders agree with" is replaced by "In the forward limit the one-loop remainders agree
  with". This relation is established numerically, not derived.
- **Sec. 4.4 and Fig. 4.** The quoted forward coefficients are determined from four points
  (cubic polynomial, M₁) and three points (quadratic polynomial, M₂); the text no longer calls
  this a fit. The comparison with one more point is described as a least-squares fit. Fig. 4
  caption: "fitted over" is replaced by "determined from".
- **Conclusions.** "including forward angles" is replaced by "including points close to the
  beam axis (|c| ≤ 0.999)", and the forward-limit agreement is stated as exact, with a reference
  to eq. (4.12).
- **App. D.** The reproducibility sentence now says that Fig. 1 is generated from its diagram
  source and that the gg → HH check starts from the supplied coefficients.
- **References.** Gröber and Mühlleitner are spelled with umlauts in the entry for
  arXiv:1705.05314.
- **Package version 2.1.0** (`CITATION.cff`, `ggaa.__version__`). After merging, publish a
  GitHub release v2.1.0 so that Zenodo archives this version under the concept DOI.

## References added on 6 November 2026

- **Sec. 1.** The papers of Blok, Dokshitzer, Frankfurt and Strikman are cited:
  Phys. Rev. D 83 (2011) 071501, Eur. Phys. J. C 72 (2012) 1963 and Eur. Phys. J. C 74 (2014)
  2926. A new sentence before the description of the work of Gaunt and Stirling states that
  they found the one-loop contribution to four-jet production in which both incident partons
  split to be enhanced only logarithmically when the transverse momenta of the two jet pairs
  balance. The three papers are also cited for the contribution in which only one of the two
  parton pairs arises from a perturbative splitting.
