"""Access to the 1080-record comparison with the unexpanded remainders."""

import csv
from fractions import Fraction

import mpmath as mp

from .normalization import FILE_FACTOR
from .paths import VALIDATION

TABLE = VALIDATION / "validation_1080.csv"


def records(path=TABLE):
    """Rows of the validation table with exact kinematics and mpmath complex values.

    The CSV stores the values in the normalization of the published 2019 files; they are
    returned divided by 64 (ggaa.normalization)."""
    with open(path) as handle:
        for row in csv.DictReader(handle):
            yield {
                "quantity": row["quantity"],
                "mA": int(row["mA_GeV"]), "muR": int(row["muR_GeV"]), "Nc": int(row["Nc"]),
                "beta": Fraction(row["beta"]), "costheta": Fraction(row["costheta"]),
                "exact": _complex(row["re_exact"], row["im_exact"]),
                "series4": _complex(row["re_series4"], row["im_series4"]),
                "series6": _complex(row["re_series6"], row["im_series6"]),
            }


def _complex(re, im):
    """The table keeps 40 significant digits; read them without rounding."""
    with mp.workdps(45):
        return mp.mpc(re, im) / FILE_FACTOR


def relative_error(approx, exact):
    return abs(approx - exact) / abs(exact)
