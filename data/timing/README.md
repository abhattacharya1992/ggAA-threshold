# Timing benchmark

Measurements of the Mathematica timing run of 2026-09-25 (host and software in
`timing_run.json`). `timing_benchmark.csv` holds, for each point, the first-call and warm
(median of `WarmRepeats`) evaluation times of the two complete two-loop coefficients, from the
unexpanded expressions and from the sixth-order series, at `WorkingPrecision` digits.
`timing_run.json` also records the on-disk sizes of the two representations that were timed.
The manuscript table and all quoted timing numbers are generated from these files
(`scripts/make_tables.py`, `scripts/make_numbers.py`).
