"""Normalization of the two-loop remainders.

The package uses one normalization at every order: the amplitude is
C_G^2 [M^(0) + a_s M^(1) + a_s^2 M^(2)] delta^{ab}, with the tree amplitude of
ggaa.born and the one-loop remainders of arXiv:1909.08993 as published.

The published two-loop remainders (ancillary files Finite_2Loop_P1, _P2) are 64 times
larger: the 2019 two-loop renormalization step multiplied the O_G-O_G amplitude,
including the tree and one-loop terms of the infrared subtraction, by 64
(data/provenance_2019/README.md).  Everything derived from those files, the
construction files in data/construction and the validation table in data/validation,
is therefore divided by FILE_FACTOR when it enters the package.
"""

FILE_FACTOR = 64
