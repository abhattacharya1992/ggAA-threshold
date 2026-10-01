# Two-loop gg → AA amplitudes near threshold and the double parton scattering singularity

This repository contains the threshold expansion of the two-loop Type II finite
remainders M₁ and M₂ for gg → AA in the heavy top quark effective theory, the data used
to validate it, and the code that produces every number, table and figure of the
accompanying paper.

The exact remainders were published with arXiv:1909.08993 [JHEP 02 (2020) 121].
Here they are expanded about the pair threshold β_A → 0 at fixed scattering angle,
through β_A⁶, for each coefficient of n_f⁰, n_f¹ and n_f². The final series in
`data/series` take 40 kB, against 204 MB for the published two-loop files.
The published two-loop files are 64 times larger than the tree and one-loop
amplitudes, and the published one-loop files omit their μ_R dependence (they are the
remainders at μ_R = 1 GeV). Everything here is in one consistent normalization:
two loops = published ÷ 64, one loop = published + 2β₀ ln(μ_R²/GeV²) M⁽⁰⁾.

## Contents

```
data/construction/   coefficients of β^n (n = -4 ... 6) as produced by the series construction
data/runs/           outputs of the construction and validation runs
data/series/         final series M_i = c0 + β² c2 + β⁴ c4 + β⁶ c6, one file per n_f power
data/validation/     1080 comparisons with the unexpanded remainders; one-loop reference values
data/constants/      the thirteen endpoint constants MB[p,q,a,w]
src/python/ggaa/     series evaluator, exact reduction, Li_{2,2}, evaluator of the published remainders
src/mathematica/     Wolfram Language evaluator of the final series
src/construction/    code of the threshold construction and of the validation (Wolfram Language, C++, Python)
scripts/             one script per result of the paper
tests/               quick consistency tests
results/             output of the scripts
manuscript/          LaTeX source; figures/ and tables/ are written by the scripts
inputs/              download script and checksums for the arXiv:1909.08993 ancillary files
```

## Requirements

Python 3.9 or later with `sympy >= 1.12`, `mpmath`, `numpy` and `matplotlib`
(`pip install -r requirements.txt`). The manuscript needs a TeX installation with
REVTeX 4.2. A Wolfram kernel is needed only for `src/mathematica` and for rerunning the
construction (`make reproduce`, which also needs g++ with GMP or Boost.Multiprecision).

## Quick start

```bash
python3 scripts/evaluate.py --beta 1/10 --costheta 1/3
```
```
M1 = -254502447.9313525400491509 - 56358840.57431776321684063 i
M2 = 6604492.650225327756854782 + 8109371.934730271156083198 i
```

Options: `--mA`, `--muR` (GeV), `--Nc`, `--nf`, `--order 4|6`, `--digits`.
From Python:

```python
from ggaa.series import ThresholdSeries        # with src/python on the path
series = ThresholdSeries()
series(1, "1/10", "1/3", mA=125, muR=125, Nc=3, nf=5, order=6)
series.component(2, 1, "1/10", "1/3")          # coefficient of n_f^1 in M_2
```

`make test` runs the tests in about ten seconds.

## Reproducing the paper

| result | command | output |
|---|---|---|
| inverse and odd powers of β vanish; final series | `make series` | `data/series/`, `results/reduce_series.txt` |
| reductions of the endpoint constants (App. A) | `make constants` | `data/constants/mb_constants.csv`, `results/mb_constants.txt` |
| leading power in terms of ℒ | `make leading-logs` | `results/leading_logs.txt` |
| one-loop threshold and forward limits | `make inputs one-loop` | `results/one_loop.txt` |
| independent evaluation of the exact two-loop remainders | `make inputs exact-check` | `results/exact_amplitude.txt` |
| Table III, error slopes | `make validation` | `manuscript/tables/accuracy.tex`, `results/compare_validation.txt` |
| Figs. 2, 3 and 6 (validation, DPS weight, flavour), Fig. 1 (diagrams); interference and ratio plots for the supplement | `make figures` | `manuscript/figures/fig1, fig7, fig2, fig0 (fig3, fig4)` |
| Figs. 4 and 5 (amplitudes by order, cross section) | `make inputs virtual physics-figures` | `manuscript/figures/fig5, fig6` |
| one-loop ln p_T coefficient against the DPS helicity weight F_A(x) | `make inputs dps` | `results/dps_check.{txt,json}` |
| ln p_T terms of V1, V2; coefficient −(7/4)π²N_c²; v0, v1 at three angles and two scales; maximum of V2 | `make inputs rate-logs` | `results/rate_logs.json` |
| velocity below which the two-loop term of M1 exceeds the one-loop term (Fig. 4) | `make crossover` | `results/amplitude_crossover.json` |
| one-loop gg → HH forward logarithm against F_S; HH points of Fig. 3 (data in `data/dps_hh`) | `make dps-hh` | `results/dps_hh.json` |
| earlier 2019 run vs published renormalisation, June 2019 pole record (App. B) | `make renorm-2019` | `results/renorm_2019.json` |
| construction files and validation records against the outputs of the construction runs | `make provenance` | `results/provenance.json` |
| inputs of the P1 n_f construction against the published remainder | `make nf1-inputs` | `results/nf1_inputs.txt` |
| construction rerun from the published files (Mathematica 13.2, several hours) | `make inputs reproduce` | `build/reproduce/` |
| diagrams of Fig. 1 | `make diagrams` | `manuscript/figures/fig0_diagrams_tikz.tex`, fallback `fig0_diagrams.pdf` |
| share of the ln p_T terms in V2 at fixed angle | `make inputs log-share` | `results/log_share.json` |
| 2019 archive against the published files (one loop: mu_R term; two loops: factor 64) | `make inputs archive` | `results/archive_2019.json` |
| top quark mass factor of the Type II tree amplitude at threshold | `make top-mass` | `results/top_mass.json` |
| comparison curves for Fig. 4 of arXiv:1909.08993 | `make fig4-2019` | `results/fig4_2019_check.{txt,pdf}` |
| forward-angle and constant tables | `make tables` | `manuscript/tables/` |
| every number quoted in the text | `make inputs numbers` | `manuscript/numbers.tex` |
| LO partonic cross section (Type II) | `make lo` | `results/lo_cross_section.{txt,csv,pdf}` |
| finite virtual corrections V1, V2 | `make inputs virtual` | `results/virtual_correction.{txt,csv,pdf}` |
| normalization of tree, one- and two-loop amplitudes | `make inputs normalization` | `results/normalization.txt` |
| paper | `make inputs paper` | `manuscript/main.pdf` and the layout report |
| arXiv submission (after `make paper`) | `make arxiv` | `build/ggAA-threshold-arxiv.tar.gz` |

`make all` needs the downloaded inputs (`make inputs`) for the numbers of the text; it runs everything
except the one-loop and exact-amplitude checks; `make exact-check` takes about 45 minutes. `make constants` takes a few minutes at 110 digits; the other targets take
at most a few minutes each.

## Conventions

- β_A = (1 − 4m_A²/ŝ)^{1/2}; cos θ is the angle between p₁ and p₃ in the partonic
  centre-of-mass frame; the kinematics are those of Sec. II B of the paper.
- Overall couplings and Wilson coefficients are stripped, as in the ancillary files of
  arXiv:1909.08993. The results are finite remainders, not cross sections.
- The series files use `beta`, `costh`, `m2` = m_A², `mu2` = μ_R², `nc` = N_c, and the
  functions `Log`, `PolyLog[2, (1 + costh)/2]`, `PolyLog[4, 1/2]`, `Zeta[3]`.
- The expansion is at fixed angle. It was compared with the unexpanded remainders for
  10⁻³ ≤ β_A ≤ 0.3 and |cos θ| ≤ 2/3 on three parameter sets, and for β_A = 0.1, 0.2, 0.3
  up to |cos θ| = 0.999 (`data/README.md`). Towards the beam the truncation error stays
  within the range found on the main grid.

## Main results

All two-loop numbers are in the common normalization (published files ÷ 64).

- After the reduction of the endpoint constants, every inverse power β⁻⁴ … β⁻¹ and
  every odd power β, β³, β⁵ vanishes identically, for arbitrary angle, m_A, μ_R and N_c,
  and every coefficient is even in cos θ.
- Define ℒ = ½ ln[−4β_A²(1 − cos²θ)] = ln(2p_T/m_A) + iπ/2 + O(β_A²). At β_A⁰ the angle
  enters only through ℒ:

  M₁ = π m_A² [ −(28/3) N_c² ℒ³ − (14/3) N_c (11N_c − 2n_f) ℒ² + … ],
  M₂ = π m_A² [ −2 N_c² ℒ² − (4/3) N_c (4N_c − n_f) ℒ + … ] + O(β_A²).

- At one loop, M₁ = m_A² (7π N_c ℒ + const) and M₂ → (π/2) N_c m_A² at threshold, both
  independent of n_f and μ_R.
- The ln p_T term is the double parton scattering singularity (Gaunt and Stirling,
  arXiv:1103.1888). In the forward limit the one-loop coefficient of ln p_T² equals
  π N_c m_A² F_A(x), with F_A(x) = 1/(x(1−x)) − 2x(1−x) and x = (1+β_A)/2, to 3×10⁻¹⁵
  for 0.001 ≤ β_A ≤ 0.7 (`make dps`). A scalar would have + instead of −.
- In the finite virtual corrections, V₁ has no ln p_T at threshold; V₂ contains
  ln² p_T with coefficient −(7/4)π²N_c² (`make rate-logs`).

## Manuscript data

No number in the manuscript is typed by hand. `make numbers` computes every quoted value from
`data/` and `results/` and writes `manuscript/numbers.tex`; `make tables`, `make validation` and
`make leading-logs` write the tables, and `make figures` and `make physics-figures` the figures.
The benchmark parameters are set once in `src/python/ggaa/benchmark.py`. `make paper` regenerates
tables and numbers before compiling.

## Physical quantities

`src/python/ggaa/born.py` gives the tree-level Type II amplitude, M_1^(0) = i ŝ/2 and
M_2^(0) = i m_A² (x − z)(x² + 1)(xz − 1)/(4xz(x² − xz + 1)), in the normalization of the loop
remainders (`data/born/README.md`). It agrees with a direct evaluation of the two Feynman
diagrams, which fixes the amplitude as C_G² M_i δ^{ab} with
C_G = −2^{5/4} a_s G_F^{1/2} cot β_H. `ggaa.crosssection` builds from it the helicity sums, the
LO partonic cross section and the finite virtual corrections

  Σ|M|² = C_G⁴ Σ_i |M_i^(0)|² [1 + a_s V₁ + a_s² V₂],

where V₁ and V₂ are Catani-subtracted finite parts, not observables.

Three points limit their interpretation:

- In the effective theory of Eq. (3) the Type II diagrams are the whole tree amplitude. In a
  model with an AAgg contact term or scalar exchange (Type I), Type I enters at lower order in
  a_s and dominates; the Type II cross section is then a small subset (≈10⁻⁸ pb for
  cot β_H = 1 near threshold). The cross section also grows with ŝ, as expected for the
  effective operators.
- V₂ uses the two-loop threshold series and is limited to β_A ≲ 0.3.
- All orders use one normalization: amplitude = C_G² [M⁽⁰⁾ + a_s M⁽¹⁾ + a_s² M⁽²⁾] δ^{ab}.
  Two corrections to the published files are needed (`data/provenance_2019/README.md`):
  (i) the published two-loop remainders are 64 times larger: the 2019 two-loop
  renormalization step multiplied the O_G–O_G amplitude, including the tree and one-loop terms
  of the infrared subtraction, by 64. The factor is removed in one place,
  `src/python/ggaa/normalization.py`: `data/series` is written divided by 64, and the
  validation table (stored in the published normalization) is divided by 64 when read.
  (ii) The published one-loop remainders do not depend on μ_R and change under a common
  rescaling of m_A and μ_R (they contain ln m_A² where ln(m_A²/μ_R²) belongs).
  `ggaa.oneloop.OneLoop` adds the missing term 2β₀ ln(μ_R²/GeV²) M⁽⁰⁾ by default
  (`scale_fix=False` returns the published expression).
  The 2019 archive in `data/provenance_2019/archive` (the original Mathematica runs) confirms
  both: its one-loop output has the mu_R term, and the published two-loop P1 remainder is 64 times
  its two-loop output (`make archive`).
  `make normalization` checks both: the corrected one-loop remainders scale correctly, and
  dM⁽²⁾/d ln μ² = 3 β₀ M⁽¹⁾ + (2β₁ − 8H_g⁽²⁾ − 12 C_F n_f) M⁽⁰⁾ holds to the accuracy of the
  series (the published files give 192 β₀ M⁽¹⁾). Relative errors, cancellations, the
  structure of the expansion and all ln p_T coefficients do not depend on either correction.
- V₂ does not include the O_G–O_J contribution proportional to the physical C_J (m_t dependent).

## Provenance

`data/construction` was produced by the code in `src/construction`, which expands the exact
remainders of arXiv:1909.08993 about threshold, and the reference values in `data/validation`
by a separate run of that code that evaluated the unexpanded remainders at 80–300 digits and
compared them with the saved series (`src/construction/README.md`). The outputs of these runs
are kept in `data/runs`; `make provenance` checks that every construction file and every
validation record agrees with them. `make reproduce` repeats the construction from the
published files and compares the result with `data/construction`. Everything downstream of
`data/construction` and `data/validation` is reproduced by this repository with Python alone.

The reference values are checked independently: `ggaa.twoloop` evaluates the published
two-loop files directly, with its own parser and an implementation of Li_{2,2} from its
defining double sum (`ggaa.li22`), and shares no code with the Wolfram Language
implementation. `make exact-check` compares the two at four points covering all three
parameter sets, and the flavour decomposition at one of them; all agree to the 40 digits
stored in the table (`results/exact_amplitude.txt`). The published expressions contain
removable 0/0 terms at cos θ = 0 (t = u), so they cannot be evaluated there directly.

## Use of AI tools

AI tools (Claude, Anthropic) were used to help write and debug the code of the
construction and validation (`src/construction`), to help write and check parts of the
analysis and plotting code, and to reorganize this repository. The author defined the
calculations, ran the code and checked the results reported in the paper.

## Licence and citation

Code: BSD 3-Clause (`LICENSE`). Data and figures: CC BY 4.0. The exact remainders of
arXiv:1909.08993 are not redistributed; `inputs/fetch_inputs.sh` downloads them. The one exception
is the coefficient of n_f in the P1 two-loop remainder, `src/construction/inputs/P1_nf1.wl` and its compact
form, which the construction reads. The code used for
the one-loop gg → HH amplitudes (`data/dps_hh`) is not included; it is available from the author on request.
Please cite the paper; see `CITATION.cff`.

## Building the paper and checking the layout

Unzip each release into an **empty** directory (unzipping over an older copy can leave stale files).
`make paper` builds `manuscript/main.pdf` with RevTeX and then runs `scripts/check_layout.py`, which
prints a plain-text report from `main.aux` and `main.log`: the page of every figure and table against
the pages of the section that first cites it, stuck floats, overfull or stretched columns, and whether
floats were flushed onto a page of their own before the references. Without the class option
`floatfix` (removed in 6.0.0) RevTeX's "A float is stuck" notices are information only: RevTeX prints
them whenever a single-column float waits at a page break, and the table shows where the float went.
With pdftotext installed the report also names columns with few words (SPARSE), which may hold only
floats and are worth a look in the PDF. The first line gives the package
version, so an old copy is recognised at once. `make layout` reruns the report without rebuilding.

`scripts/check_layout_standin.py` is a stress test for sandboxes without RevTeX: it builds an
article-class look-alike with 12 shifts of all page breaks and applies the same report. It is not a
preview of the RevTeX layout.

## Figures and tables

All plots are drawn by the scripts from the data files, in the usual journal style: boxed axes with
inward major and minor ticks on all sides, Times-type labels matching the text, black and standard
colours always paired with distinct line styles or open symbols, frameless legends inside the frame
and no titles. Single-column figures are 3.375 in wide and two-column figures 7.0 in, with lettering
of at least 7.5 pt. Figure 1 is drawn with TikZ-Feynman from one source, `manuscript/figures/fig0_diagrams_tikz.tex`
(written by `scripts/make_diagrams.py`). When `tikz-feynman.sty` is installed (`kpsewhich tikz-feynman.sty`),
the manuscript sets the diagrams directly, in the fonts of the paper (`manuscript/diagram_setup.tex`);
otherwise it includes the prebuilt `fig0_diagrams.pdf` compiled from the same source. Tables use RevTeX's `ruledtabular`.

## Submitting to arXiv

After `make paper`, `make arxiv` writes `build/ggAA-threshold-arxiv.tar.gz`: the manuscript sources,
figures and `main.bbl` (arXiv does not run BibTeX), and this package in `anc/` as ancillary files.
Following arXiv's rules `anc/` contains no TeX files, so it holds everything except `manuscript/`;
the downloaded files of arXiv:1909.08993 are left out (`inputs/fetch_inputs.sh` fetches them).
