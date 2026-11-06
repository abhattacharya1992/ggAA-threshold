(* Evaluate the threshold series of data/series in Wolfram Language.

   In a notebook or script:
       Get["src/mathematica/evaluate_series.wl"];
       ThresholdM[1, 1/10, 1/3]                      (* M_1 at beta = 1/10, cos(theta) = 1/3 *)
       ThresholdM[2, 1/10, 1/3, "muR" -> 250, "Order" -> 4]

   From the shell, `wolframscript -file src/mathematica/evaluate_series.wl` prints the
   benchmark point.  Expected output, through beta^6 with m_A = mu_R = 125 GeV, N_c = 3, n_f = 5:
       M1 = -2.545024479313525400491509 10^8 - 5.635884057431776321684063 10^7 I
       M2 =  6.604492650225327756854782 10^6 + 8.109371934730271156083198 10^6 I
   (the series files hold the published two-loop remainders divided by 64; see README).   *)

seriesDirectory = FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "data", "series"}];

seriesCoefficients[i_Integer, k_Integer] := seriesCoefficients[i, k] =
  Get[FileNameJoin[{seriesDirectory, "M" <> ToString[i] <> "_nf" <> ToString[k] <> ".m"}]];

Options[ThresholdM] = {"mA" -> 125, "muR" -> 125, "Nc" -> 3, "nf" -> 5, "Order" -> 6, "Digits" -> 30};

ThresholdM[i : 1 | 2, b_, c_, OptionsPattern[]] := Module[{n = OptionValue["nf"], powers, total},
  powers = Range[0, OptionValue["Order"], 2];
  total = Sum[n^k Sum[beta^p seriesCoefficients[i, k][[p/2 + 1]], {p, powers}], {k, 0, 2}];
  N[total /. {beta -> b, costh -> c, m2 -> OptionValue["mA"]^2, mu2 -> OptionValue["muR"]^2,
              nc -> OptionValue["Nc"]}, OptionValue["Digits"]]
];

If[Length[$ScriptCommandLine] > 0,
  Print["M1 = ", ThresholdM[1, 1/10, 1/3]];
  Print["M2 = ", ThresholdM[2, 1/10, 1/3]]
];
