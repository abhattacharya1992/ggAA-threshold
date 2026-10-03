# Construction of the threshold series

Code that expands the published two-loop Type II remainders of arXiv:1909.08993 about
threshold and evaluates the validation records. It was run in September 2026 with
Mathematica 13.2; the outputs of those runs are in `data/runs`.

| folder | content |
|---|---|
| `nf2/` | coefficients of n_f²: extraction from the published files (`extract_nf.cpp`, `check_nf_extraction.cpp`), compact form and checks (`parse_nf2.py`, `collect_nf2_normalized.py`, `derive_nf2.py`, `validate_nf2_python.py`), expansion (`driver.wl`) |
| `p1_nf1/construction/` | coefficient of n_f in M1: expansion of the 526 blocks of `blocks.wl` through β⁴ and β⁶ |
| `p1_nf1/check/` | exact checks of the saved P1 n_f series against an independent recurrence (`python_reference_series.wl`) and numerical checks at 135 points |
| `remaining/` | coefficients of n_f⁰ in M1 and M2 and of n_f in M2: collection of the published expressions into blocks (`first_run/collect_complete.cpp`, `first_run/check_collection.cpp`), exact rational series of the blocks in β (`prepare_rational_series.py` with `rational_dag.cpp`, `jet_dag.cpp`, `reference_dag.cpp`, `check_fixed_reference.cpp`, `exact_rf.hpp`), expansion of the special functions (`UnifiedMPL.wl`) and assembly of M1 and M2 (`driver.wl`) |
| `validation/` | evaluation of the unexpanded remainders and comparison with the saved series at 1080 records (`validate_saved.wl`) |
| `inputs/` | the coefficient of n_f in the published P1 remainder (`P1_nf1.wl`) and its compact form (`P1_nf1_compact.wl`) |

The special functions are expanded as in Sec. III A of the paper: their derivatives in β are
expanded and integrated term by term from their values at β = 0. The endpoint constants
`MB[p,q,a,w]` are the integrals of Eq. (A1); their reduction to the basis of eight constants
and the exact cancellations are done by `make constants` and `make series`.

In September the remaining coefficients were built in two runs: the first prepared the
rational series of P1 blocks 0 to 1437 in Mathematica and stopped at block 1438, and the second
prepared the other blocks with the C++ programs and completed the construction.
`reproduce_series.sh` prepares all blocks with the C++ programs. The code here is the code of
those runs with its comments, formatting, variable names and printed messages edited; the
calculation is unchanged.

## Running

From the package root, with the published files downloaded (`make inputs`):

    make reproduce

which runs

    bash src/construction/reproduce_series.sh --inputs inputs/arXiv-1909.08993/anc \
        --nf1-inputs src/construction/inputs --work build/reproduce --package .

It needs Mathematica 13.2 (`math`, `WolframKernel` or `wolfram` on the path), g++ with C++17
and GMP or Boost.Multiprecision, and Python 3.9 or later with sympy and mpmath. It checks the
SHA-256 of all four input files first, skips stages that have finished when it is run again,
and ends by comparing the result with every file in `data/construction`. The construction
driver of the P1 n_f coefficient stops after writing its two series, as it did in September;
the following check stage then verifies them. On lhcpheno48, with Boost.Multiprecision, the whole run took about 4.5 hours with a peak
memory of 0.8 GB; the exact rational series of the 14,247 blocks take most of that time. If the
GMP and Boost headers are not installed, unpack Boost (header-only) in the home directory and set
`CPLUS_INCLUDE_PATH` to it before starting. `scripts/watch_reproduce.sh` starts the run detached
from the terminal and stops it if memory runs short. The validation (about ten hours) is not
repeated by the script.

`make nf1-inputs` checks that `inputs/P1_nf1.wl` is the coefficient of n_f in the published P1
remainder and that the compact form agrees with it.
