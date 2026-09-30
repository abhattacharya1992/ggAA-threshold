# Wolfram Language

`evaluate_series.wl` evaluates the final series in `data/series`:

```
wolframscript -file src/mathematica/evaluate_series.wl
```

prints M₁ and M₂ at β_A = 1/10, cos θ = 1/3, which must agree with
`python3 scripts/evaluate.py --beta 1/10 --costheta 1/3`. Inside a session,
`Get` the file and call `ThresholdM[i, beta, costheta, "mA" -> 125, "muR" -> 125,
"Nc" -> 3, "nf" -> 5, "Order" -> 6]`.
