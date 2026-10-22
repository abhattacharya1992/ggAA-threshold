# Two-loop gg → AA amplitudes near threshold and the double parton scattering singularity

Paper: A. Bhattacharya, *Two-loop gg → AA amplitudes near threshold and the double parton
scattering singularity*, [arXiv:2610.19024](https://arxiv.org/abs/2610.19024) [hep-ph].

```bibtex
@article{Bhattacharya:2026ggAA,
    author        = "Bhattacharya, Arunima",
    title         = "{Two-loop $gg\to AA$ amplitudes near threshold and the double parton scattering singularity}",
    eprint        = "2610.19024",
    archivePrefix = "arXiv",
    primaryClass  = "hep-ph",
    year          = "2026"
}
```

Versions: the tag `arxiv-v1` (commit 097cce7) is the state submitted as arXiv v1; its
`manuscript/main.pdf` is the arXiv v1 PDF. The `main` branch is the current state. The older tags
v7.2.0 to v9.0.0 are development snapshots from before the submission.

This repository contains the threshold expansion of the two-loop Type II finite
remainders M₁ and M₂ for gg → AA in the heavy top quark effective theory, the data
used to validate it, and scripts for reproducing the numerical results, tables
and figures of the accompanying paper.

The exact remainders were published with
[arXiv:1909.08993](https://arxiv.org/abs/1909.08993), JHEP 02 (2020) 121.
Here they are expanded about the pair threshold β_A → 0 at fixed scattering
angle through β_A⁶, separately for the coefficients of n_f⁰, n_f¹ and n_f².
The final series in `data/series` occupy 40 kB, compared with 204 MB for the
published two-loop files.

The released two-loop remainders carry an extra overall factor 64 relative to
the normalization used here. The released one-loop files omit their μ_R
dependence and correspond to μ_R = 1 GeV. The package puts both orders in a
common normalization:

```text
two loops = published ÷ 64
one loop  = published + 2β₀ ln(μ_R²/GeV²) M⁽⁰⁾
β₀        = (11N_c − 2n_f)/3
```

The recommended range of the series is 0 < β_A ≤ 0.2 at the validated angles.
At the tested points in this range, the maximum relative errors of the two
complete remainders are 3.23×10⁻⁶ and 6.77×10⁻⁵, respectively. These are
tested-point results, not a uniform error bound over all kinematics.

## Contents

```text
data/construction/   coefficients of β^n, n = −4 ... 6, from the series construction
data/runs/           outputs of the construction and validation runs
data/series/         final series c0 + β² c2 + β⁴ c4 + β⁶ c6, one file per n_f power
data/validation/     1080 comparisons with the unexpanded remainders; one-loop reference values
data/constants/      numerical values and reductions of the thirteen endpoint constants
data/born/           tree-level amplitudes and checks of their normalization
data/dps_hh/         supplied one-loop gg → HH coefficients used in the DPS comparison
data/provenance_2019/ original calculation records used for the normalization checks
src/python/ggaa/     series evaluation, symbolic reduction and evaluation of the published remainders
src/mathematica/     Wolfram Language evaluation of the final series
src/construction/    threshold construction and validation code in Wolfram Language, C++ and Python
scripts/             scripts for the numerical results, tables, figures and checks
tests/               quick consistency tests
results/             saved outputs of the scripts
manuscript/          LaTeX source, tables and figures
inputs/              download script and checksums for the arXiv:1909.08993 ancillary files
```

The coefficients c0, c2, c4 and c6 include logarithms of β_A. The complete
two-loop remainder is the sum of the three flavour coefficients weighted by
1, n_f and n_f².

## Requirements

Python 3.9 or later with SymPy, mpmath, NumPy and Matplotlib. Install the
dependencies from the package root:

```bash
python3 -m pip install -r requirements.txt
```

The manuscript uses the standard LaTeX `article` class with the JHEP style.
Building it requires PDFLaTeX, BibTeX and `latexmk`. The files `jheppub.sty`
and `JHEP.bst` are supplied in `manuscript/`. TikZ-Feynman is optional because
a compiled diagram PDF is supplied as a fallback.

A Wolfram kernel is needed for the Wolfram Language implementation and for
rerunning the construction. The construction was tested with Mathematica
13.2 and also requires g++ with C++17 and either the GMP development headers
and libraries or the Boost.Multiprecision headers. The construction script
uses GMP when available and otherwise uses Boost. Boost is header-only; a
local installation is sufficient if its root is on `CPLUS_INCLUDE_PATH`.
The construction and memory-monitoring scripts use Bash and Linux utilities.

## Quick start

Run commands from the package root. Evaluating the supplied series does not
require the downloaded exact remainders:

```bash
python3 scripts/evaluate.py --beta 1/10 --costheta 1/3
```

```text
M1 = -254502447.9313525400491509 - 56358840.57431776321684063 i
M2 = 6604492.650225327756854782 + 8109371.934730271156083198 i
```

Options are `--mA`, `--muR` in GeV, `--Nc`, `--nf`, `--order 4|6` and
`--digits`. The defaults are m_A = μ_R = 125 GeV, N_c = 3, n_f = 5,
order 6 and 30-digit working precision.

From Python:

```python
import sys
import mpmath as mp

sys.path.insert(0, "src/python")
from ggaa.series import ThresholdSeries

mp.mp.dps = 50
series = ThresholdSeries()
value = series(1, "1/10", "1/3", mA=125, muR=125, Nc=3, nf=5, order=6)
component = series.component(2, 1, "1/10", "1/3")  # coefficient of n_f¹ in M₂
```

The Wolfram Language implementation evaluates the same supplied series:

```bash
wolframscript -file src/mathematica/evaluate_series.wl
```

`make test` runs the quick tests. Tests that need the published one-loop
files are skipped when those files are absent. To download the inputs and
run the complete suite:

```bash
make inputs test
```

## Reproducing the paper

The commands below use supplied data and saved results where indicated.
`make inputs` downloads the original ancillary files and verifies their
checksums. Targets that need these files include it explicitly in the table.

| Result or check | Command | Output |
|---|---|---|
| Inverse and odd powers cancel; final series | `make series` | `data/series/`, `results/reduce_series.txt` |
| Endpoint-constant reductions, App. A | `make constants` | `data/constants/mb_constants.csv`, `results/mb_constants.txt` |
| Leading-power coefficients in terms of ℒ | `make leading-logs` | `results/leading_logs.txt`, `manuscript/tables/leading_logs.tex` |
| One-loop threshold and forward limits | `make inputs one-loop` | `results/one_loop.txt` |
| Independent evaluation of the exact two-loop remainders | `make inputs exact-check` | `results/exact_amplitude.txt` |
| Table 3 and error slopes | `make validation` | `manuscript/tables/accuracy.tex`, `results/compare_validation.txt` |
| Figs. 1, 2, 4 and 6; additional interference and ratio plots | `make figures` | `manuscript/figures/` |
| Figs. 3 and 5, amplitudes by order and Type II cross section | `make inputs virtual physics-figures` | `manuscript/figures/fig5_amplitudes.pdf`, `manuscript/figures/fig6_cross_section.pdf` |
| One-loop DPS logarithm against the helicity weight F_A(x) | `make inputs dps` | `results/dps_check.txt`, `results/dps_check.json` |
| Logarithms of V₁ and V₂, including the coefficient −(7/4)π²N_c² | `make inputs rate-logs` | `results/rate_logs.json` |
| Velocity at which the two-loop term of M₁ exceeds the one-loop term, Fig. 3 | `make inputs crossover` | `results/amplitude_crossover.json` |
| Supplied one-loop gg → HH forward coefficients against F_S, Fig. 4 | `make dps-hh` | `results/dps_hh.json` |
| Earlier 2019 renormalisation and the June 2019 pole record, App. B | `make renorm-2019` | `results/renorm_2019.json` |
| Opposite-helicity forward fits supporting the real ℓ² coefficient −4N_c times the one-loop forward remainder, Sec. 4.4 | `make inputs jz2` | `results/jz2_forward_nc3_nf5_e567.txt`, `results/jz2_forward_nc4_nf2_e567.txt` |
| Two-loop logarithmic coefficient relations through β_A⁶, including terms beyond leading power, Sec. 4.4 | `make nlp-dps` | `results/nlp_dps.txt` |
| Equal-helicity forward fits at the tested velocities in 0.05 ≤ β_A ≤ 0.7, Sec. 4.4 and Fig. 4 | `make inputs m1-forward` | `results/m1_forward_Nc3_nf5.txt`, `results/m1_forward_Nc4_nf2.txt` |
| Construction files and validation records against the archived run outputs | `make provenance` | `results/provenance.json` |
| Supplied P1 n_f construction inputs against reference values and each other | `make nf1-inputs` | `results/nf1_inputs.txt` |
| Rerun the threshold construction | `make inputs reproduce` | `build/reproduce/`, `build/reproduce/reproduce.log` |
| Diagrams of Fig. 1 | `make diagrams` | `manuscript/figures/fig0_diagrams_tikz.tex`, `manuscript/figures/fig0_diagrams.pdf` |
| Share of the logarithmic terms in V₂ at fixed angle | `make inputs log-share` | `results/log_share.json` |
| Original archive against the published files, including the scale term, factor 64 and counterterm comparison | `make inputs archive` | `results/archive_2019.json` |
| Top quark mass factor of the Type II tree amplitude at threshold | `make top-mass` | `results/top_mass.json` |
| Comparison curves for Fig. 4 of arXiv:1909.08993 | `make fig4-2019` | `results/fig4_2019_check.txt`, `results/fig4_2019_check.pdf` |
| Forward-angle and endpoint-constant tables | `make tables` | `manuscript/tables/` |
| Numerical macros used in the manuscript | `make inputs numbers` | `manuscript/numbers.tex` |
| LO Type II partonic cross section | `make lo` | `results/lo_cross_section.txt`, `.csv`, `.pdf` |
| Tree-weighted finite virtual corrections | `make inputs virtual` | `results/virtual_correction.txt`, `.csv`, `.pdf` |
| Tree, one-loop and two-loop normalization checks | `make inputs normalization` | `results/normalization.txt` |
| Build the JHEP manuscript and print the layout report | `make inputs paper` | `manuscript/main.pdf` and the layout report |
| Prepare the arXiv submission after building the paper | `make arxiv` | `build/ggAA-threshold-arxiv.tar.gz` |

`make all` runs `series`, `constants`, `leading-logs`, `validation`, `tables`,
`figures` and `numbers`. Download the inputs first with `make inputs`.
The other targets are run separately. `make inputs all-figures` regenerates
the manuscript figures, including the virtual corrections used in Fig. 5.

The forward-fit targets normally reuse the supplied exact-amplitude points
in `results/jz2_forward_points.jsonl` and `results/m1_forward_points.jsonl`.
Missing points are evaluated from the original expressions. The equal-helicity
fits use 1 − cos θ = 10⁻⁷, …, 10⁻¹⁰; the additional `_e6-10` outputs include
10⁻⁶. Full forward scans take several hours. The one-loop gg → HH check uses
the supplied coefficients in `data/dps_hh`; it does not regenerate those
amplitudes.

`make archive` uses the saved published two-loop P1 value by default.
To recompute that value after downloading the inputs:

```bash
python3 scripts/check_archive_2019.py --compute
```

Most targets that use supplied data take seconds to a few minutes. The endpoint
constants are recomputed at 110 digits. A full uncached `make exact-check` takes
about 45 minutes on the reference machine. It saves intermediate evaluations
so that an interrupted run resumes from completed calculations. Runtimes depend
on the machine and on existing cached results.

The full construction took about 4.5 hours with 0.8 GB peak memory on
the reference machine using Mathematica 13.2 and Boost.Multiprecision. `make reproduce`
reruns the construction and compares its coefficients with the supplied files;
it does not repeat the original long validation scan. For a detached run with
memory monitoring, first download the inputs and then use, for example:

```bash
setsid nohup bash scripts/watch_reproduce.sh 60 15 > /dev/null 2>&1 < /dev/null &
```

The arguments set the maximum job memory and minimum available system memory
in GiB. Set them for the machine being used. The monitor writes `reproduce.log`
and `reproduce_mem.log` in the package root and stops the run when either limit
is reached. Finished construction stages are reused when the run is restarted.

## Conventions and validation range

- β_A = (1 − 4m_A²/ŝ)^{1/2}. The scattering angle θ is measured between p₁ and p₃
  in the partonic centre-of-mass frame, as in Sec. 2.2 of the paper.
- Overall couplings and Wilson coefficients are stripped from M₁ and M₂.
  These quantities are finite remainders.
- The series files use `beta`, `costh`, `m2` = m_A², `mu2` = μ_R² and
  `nc` = N_c, with `Log`, `PolyLog[2, (1 + costh)/2]`, `PolyLog[4, 1/2]`
  and `Zeta[3]`.
- The fixed-angle series were compared with the unexpanded remainders for
  10⁻³ ≤ β_A ≤ 0.3 and |cos θ| ≤ 2/3 on three parameter sets. Additional
  tests at β_A = 0.1, 0.2 and 0.3 extend to |cos θ| = 0.999. See
  `data/README.md` and Sec. 3.3 of the paper.
- Towards the beam, the complete-remainder errors stay within the range found
  on the main grid at the same velocity, except for M₁ at β_A = 0.1. There
  the error at cos θ = 0.999 is 21% larger than the main-grid maximum.
- The recommended range is 0 < β_A ≤ 0.2 at the tested angles. At β_A = 0.3,
  the largest relative error of a complete remainder on the validation grid
  is 1.94×10⁻³. Some plots include this point to show the behaviour beyond
  the recommended range. Use the unexpanded remainders outside the range
  recommended in the paper.

## Main results

All two-loop coefficients below use the common normalization, namely the
published files divided by 64.

- After substituting the numerically identified endpoint-constant reductions,
  every inverse power β_A⁻⁴, …, β_A⁻¹ and every odd power β_A, β_A³ and β_A⁵
  cancels symbolically. Every retained coefficient is even in cos θ.
- At leading power, the angular dependence enters only through
  ℒ = ½ ln[−4β_A²(1 − cos²θ)] = ln(2p_T/m_A) + iπ/2 + O(β_A²).
  The highest logarithmic terms are:

  ```text
  M₁⁽²⁾ = π m_A² [−(28/3)N_c² ℒ³ − (14/3)N_c(11N_c − 2n_f) ℒ² + …]
  M₂⁽²⁾ = π m_A² [−2N_c² ℒ² − (4/3)N_c(4N_c − n_f) ℒ + …]
  ```

- At one loop and leading power,
  M₁⁽¹⁾ = m_A²[7πN_c ℒ + constant] and M₂⁽¹⁾ → (π/2)N_c m_A².
  The logarithmic coefficient and the M₂ limit are independent of n_f and
  μ_R; the M₁ constant is not. The opposite-helicity forward limit is compared
  with the closed form πN_c[x⁴ + (1 − x)⁴]/[x(1 − x)] for M₂⁽¹⁾/m_A²,
  where x = (1 + β_A)/2, at the tested velocities in 0.05 ≤ β_A ≤ 0.9.
- The one-loop logarithm is the double parton scattering singularity studied
  by Gaunt and Stirling, [arXiv:1103.1888](https://arxiv.org/abs/1103.1888).
  With the overall constant fixed at threshold, the forward coefficient of
  ln(p_T²/m_A²) agrees with πN_c m_A² F_A(x), where
  F_A(x) = 1/[x(1 − x)] − 2x(1 − x), to a relative accuracy of 3×10⁻¹⁵
  at the tested velocities in 0.001 ≤ β_A ≤ 0.7 (`make dps`). The supplied
  gg → HH coefficients support the corresponding CP-even weight
  F_S(x) = 1/[x(1 − x)] + 2x(1 − x) with the same overall constant.
- The two-loop forward fits support the coefficient relations described in
  Sec. 4.4. The opposite-helicity remainder develops a double logarithm even
  though it has no one-loop DPS logarithm. These fits do not establish
  exponentiation in QCD.
- At leading power, the one-loop DPS logarithm drops out of the tree
  interference because of its phase. At the next order, combining the
  tree–two-loop interference with the one-loop square gives
  V₂ ⊃ −(7/4)π²N_c² ln²(2p_T/m_A) (`make rate-logs`).

## Manuscript data

`make numbers` generates the numerical macros in `manuscript/numbers.tex`
from the supplied data, saved results, benchmark parameters and downloaded
inputs. `make tables`, `make validation` and `make leading-logs` write the
tables. `make figures` and `make physics-figures` write the manuscript plots.
The benchmark parameters are defined in `src/python/ggaa/benchmark.py`.

After table generation, `scripts/format_jhep_tables.py` converts the tables
to the JHEP format. The manuscript Makefile runs this step before compilation.
The formatter must be present when rebuilding the JHEP manuscript from
regenerated tables.

`make paper` regenerates tables and numerical macros before compiling. It uses
the supplied figures and saved results; it does not rerun the full construction,
all validation calculations or all figure scripts. Run the relevant targets
above when regenerating those results.

## Physical quantities

`src/python/ggaa/born.py` gives the tree-level Type II tensor coefficients:

```text
M₁⁽⁰⁾ = i ŝ/2
M₂⁽⁰⁾ = i m_A² (x − z)(x² + 1)(xz − 1)/[4xz(x² − xz + 1)]
```

Here x and z are the scattering variables of Sec. 2.2. This x differs from
the momentum fraction x = (1 + β_A)/2 used in the DPS formulae above.
The tree normalization is checked against a direct evaluation of the two
Feynman diagrams (`data/born/README.md`). The helicity amplitudes have the
overall factor C_G² δ^{ab}, where
C_G = −2^{5/4} a_s G_F^{1/2} cot β_H and a_s = α_s/(4π).

`ggaa.crosssection` provides the LO Type II partonic cross section and the
finite virtual corrections. With M⁽¹⁾ and M⁽²⁾ denoting finite remainders,
the helicity and colour sum through two-loop order is:

```text
Σ_{hel,col}|A_fin|² =
2(N_c² − 1) C_G⁴ Σ_i |M_i^(0)|² [1 + a_s V₁ + a_s² V₂]
```

V₁ and V₂ are Catani-subtracted finite virtual coefficients. Real emission
and initial-state collinear factorisation are needed for a complete
higher-order cross section. The finite virtual coefficients retain subtraction
scheme dependence and are not observables or cross-section K factors.

The following points specify their use:

- The Type II tree amplitude is leading when the pseudoscalar couples to
  gluons only through the operator O_G of eqs. (1.1)–(1.2). In a 2HDM, Type I enters at a
  lower order in a_s, and a complete prediction also includes Type I and its
  interference with Type II. The Type II tree cross section is about 10⁻⁸ pb
  near threshold for the paper's benchmark parameters with cot β_H = 1.
- V₂ uses the two-loop series through β_A⁶. Its recommended input range is
  0 < β_A ≤ 0.2 at the validated angles; the plots also include β_A = 0.3.
- The amplitude is normalized as
  C_G²[M⁽⁰⁾ + a_s M⁽¹⁾ + a_s² M⁽²⁾]δ^{ab} at every order.
- V₂ does not include the mixed O_G–O_J contribution proportional to the
  physical C_J, which depends on m_t. Its treatment is described in Sec. 5.2.

The two normalization conversions are implemented as follows:

1. The original two-loop renormalisation step multiplied the O_G–O_G
   amplitude, including the tree and one-loop terms used in infrared
   subtraction, by 64. `src/python/ggaa/normalization.py` defines this factor.
   The series are written divided by 64, and the validation records, stored
   in the published normalization, are divided by 64 when read.
2. `ggaa.oneloop.OneLoop` adds the missing one-loop term
   2β₀ ln(μ_R²/GeV²) M⁽⁰⁾ by default. `scale_fix=False` returns the
   published expression. The corrected remainder has the scale dependence
   required by the renormalisation-group relation.

The records in `data/provenance_2019/archive` support both conversions.
Their one-loop outputs contain the missing μ_R dependence. The two-loop
comparison includes the factor 64 and the counterterm and residual term
discussed in App. B; it is not a simple equality between the published
remainder and 64 times the earlier output.

`make normalization` checks the one-loop scale dependence and the two-loop
relation:

```text
dM⁽²⁾/d ln μ_R² =
3β₀ M⁽¹⁾ + [2β₁ − 8H_g⁽²⁾ − 12C_F n_f] M⁽⁰⁾
```

The subtraction convention for H_g⁽²⁾ is specified in App. B. The two-loop
relation holds to the accuracy of the truncated series. Dividing the two-loop
expressions by 64 leaves their relative truncation errors, symbolic
cancellations and logarithmic powers unchanged. Restoring the one-loop scale
dependence changes its constant and V₁, V₂, but not the coefficient of the
one-loop DPS logarithm.

## Provenance

The code in `src/construction` produced the coefficients in
`data/construction`. Separate validation runs evaluated the unexpanded
remainders at 80–300 digits and compared them with the saved series.
The outputs are retained in `data/runs`. `make provenance` checks the
construction files and validation records against these outputs.

`make reproduce` reruns the construction using the downloaded remainders
and the supplied P1 n_f inputs, then compares the resulting coefficients
with `data/construction`. The construction procedure and its stages are
described in `src/construction/README.md`.

The supplied construction coefficients and validation records are sufficient
to repeat the downstream series reduction and numerical analysis in Python.
Typesetting the manuscript and compiling the TikZ diagrams require TeX.
The generation code for the separate one-loop gg → HH calculation is not
included; its coefficients are supplied for the checks used in the paper.

An independent Python implementation, `ggaa.twoloop`, evaluates the published
two-loop expressions with its own parser and an implementation of Li_{2,2}
from its defining double sum (`ggaa.li22`). It shares no implementation code
with the Wolfram Language evaluation used for the reference values.
`make exact-check` compares both complete remainders at four points covering
all three parameter sets and checks the flavour decomposition at one point.
The saved comparisons agree to the 40 digits stored in the reference table
(`results/exact_amplitude.txt`).

The published two-loop expressions contain removable 0/0 terms at
cos θ = 0, where t = u. Their direct evaluation there requires taking a
limit. The independent exact checks therefore use nonzero cos θ; the
threshold series are regular at cos θ = 0.

## Use of AI tools

The author developed and directed the calculations, analysis and software.
Large language models were used as supporting tools in the preparation,
debugging and editing of code and text, and in repository organisation.
The author is responsible for the scientific results and the released material.

The package includes comparisons with the unexpanded remainders,
renormalisation-group checks and an independent Python evaluation.
The accompanying manuscript describes the supporting use of AI tools and
the validation of the scientific calculation.

## Licence and citation

Code is released under the BSD 3-Clause licence (`LICENSE`). Data, tables and figures
are released under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).

The downloaded ancillary files of arXiv:1909.08993 are excluded from the
arXiv submission bundle. `inputs/fetch_inputs.sh` downloads them and checks
their checksums. The package includes the coefficient of n_f in the P1 two-loop
remainder, `src/construction/inputs/P1_nf1.wl`, and its compact form, which
are used by the construction. It also includes the original records needed
for the provenance checks.

The supplied one-loop gg → HH coefficients are in `data/dps_hh`. Their
generation code is not included and is available from the author on request.
Please cite the accompanying paper; see `CITATION.cff`.

## Building the paper and checking the layout

The manuscript source is `manuscript/main.tex`.

After editing the manuscript, run these commands from the package root:

```bash
(cd manuscript && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex)
```

This updates `manuscript/main.pdf`. `latexmk` runs BibTeX and additional
LaTeX passes when needed. For automatic recompilation whenever a source
file is saved:

```bash
(cd manuscript && latexmk -pdf -pvc -interaction=nonstopmode -halt-on-error main.tex)
```

Press Ctrl+C to stop watching.

To regenerate the tables and numerical macros before compilation, run:

```bash
make inputs paper
```

This also runs `scripts/check_layout.py`, which reads `main.aux` and
`main.log` and lists the pages of figures and tables relative to their
first citations. `make layout` reruns the report without rebuilding.
Its column checks were designed for the earlier two-column REVTeX layout
and do not assess the single-column JHEP layout reliably. Inspect the
compiled PDF directly, including figure and table placement.

## Figures and tables

The manuscript plots are generated from the supplied data. They use boxed
axes with inward ticks, serif labels and distinct line styles or markers.
The figure scripts set the font sizes; the LaTeX inclusion commands set
the widths in the compiled manuscript.

Figure 1 is drawn with TikZ-Feynman from
`manuscript/figures/fig0_diagrams_tikz.tex`, generated by
`scripts/make_diagrams.py`. When `tikz-feynman.sty` is installed, the
manuscript typesets these diagrams directly; otherwise it includes the
supplied `fig0_diagrams.pdf`. The JHEP tables use `tabular` with rules
from `booktabs`, with captions below the tables. The formatter described
above removes the REVTeX wrappers emitted by the table-generation scripts.

## Submitting to arXiv

The JHEP packaging script must include `jheppub.sty`, `JHEP.bst` and
`references.bib`, in addition to the manuscript sources, figures and
`main.bbl`. In `scripts/make_arxiv.py`, the initial file list inside
`manuscript_files()` is:

```python
    files, todo = {"main.tex", "main.bbl", "references.bib", "jheppub.sty", "JHEP.bst"}, ["main.tex"]
```

Compile successfully before creating the archive. From the package root:

```bash
(cd manuscript && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex) && make arxiv
```

`make arxiv` packages the files; it does not compile the manuscript or
refresh the bibliography. It writes `build/ggAA-threshold-arxiv.tar.gz`
with the manuscript sources and required figures at the archive root,
including the JHEP style files, `references.bib` and the compiled
bibliography `main.bbl`. The accompanying code and data are placed in `anc/`.
arXiv supports bibliography processing and uses a supplied `.bbl` when present.

To check the essential files:

```bash
tar -tzf build/ggAA-threshold-arxiv.tar.gz |
  grep -E '^(main\.tex|main\.bbl|references\.bib|jheppub\.sty|JHEP\.bst)$'
```

All five filenames must appear. The ancillary directory contains no TeX
files. The downloaded inputs, temporary build directories and local
evaluation caches are excluded.

Upload the submission tarball, select `main.tex` and PDFLaTeX, and check
the PDF generated by arXiv before completing the submission.

To recompile an extracted submission archive locally, run
`latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex` from its
top directory, above `anc/`. Regenerate tables and numerical macros in the
repository layout before creating a new archive so that the JHEP formatting
step is applied.

Official guidance is available for
[TeX submissions](https://info.arxiv.org/help/submit_tex.html) and
[ancillary files](https://info.arxiv.org/help/ancillary_files.html).
