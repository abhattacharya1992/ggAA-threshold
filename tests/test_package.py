"""Quick checks of the package.  Run with:  python3 -m unittest discover -s tests -v"""

import pathlib
import random
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

import mpmath as mp
import sympy as sp

from ggaa import constants, wolfram
from ggaa.paths import CONSTRUCTION, INPUTS
from ggaa.series import ThresholdSeries
from ggaa.validation import records, relative_error

SERIES = ThresholdSeries()


class EndpointConstants(unittest.TestCase):
    def test_reductions(self):
        with mp.workdps(40):
            for key, value in constants.REDUCTIONS.items():
                numeric = constants.mb_integral(*key)
                exact = mp.mpf(str(sp.N(value, 45)))
                self.assertLess(abs(numeric - exact), mp.mpf(10) ** -35, key)


class Cancellations(unittest.TestCase):
    def test_selected_inverse_and_odd_powers(self):
        for name in ("M1_nf0_beta_m2", "M2_nf0_beta_m4", "M2_nf1_beta_m1", "M1_nf0_beta_1", "M2_nf1_beta_3"):
            reduced = wolfram.reduce(wolfram.read(CONSTRUCTION / f"{name}.wl"))
            self.assertTrue(wolfram.is_zero(reduced), name)


class SeriesEvaluation(unittest.TestCase):
    def test_reference_point(self):
        # beta = 1/10, cos(theta) = 1/3, m_A = mu_R = 125 GeV, N_c = 3, n_f = 5, through beta^6.
        # The expected values are those of the published-file normalization, divided by 64.
        from ggaa.normalization import FILE_FACTOR
        with mp.workdps(40):
            m1 = SERIES(1, "1/10", "1/3")
            m2 = SERIES(2, "1/10", "1/3")
            expected1 = mp.mpc("-1.628815666760656256314565451270863831639e10",
                               "-3.606965796756336845877800095501192171822e9") / FILE_FACTOR
            expected2 = mp.mpc("4.226875296144209764387060320881623999076e8",
                               "5.189998038227373539893246631865806668357e8") / FILE_FACTOR
            self.assertLess(relative_error(m1, expected1), 1e-30)
            self.assertLess(relative_error(m2, expected2), 1e-30)

    def test_against_validation_table(self):
        rows = list(records())
        random.seed(1)
        with mp.workdps(45):
            for rec in random.sample(rows, 40):
                i, nf = int(rec["quantity"][1]), rec["quantity"].split("_nf")[1]
                b, c = rec["beta"], rec["costheta"]
                args = dict(beta=mp.mpf(b.numerator) / b.denominator, costh=mp.mpf(c.numerator) / c.denominator,
                            mA=rec["mA"], muR=rec["muR"], Nc=rec["Nc"])
                for order in (4, 6):
                    value = (SERIES(i, nf=5, order=order, **args) if nf == "5"
                             else SERIES.component(i, int(nf), order=order, **args))
                    self.assertLess(relative_error(value, rec[f"series{order}"]), 1e-35, rec["quantity"])


class TreeLevel(unittest.TestCase):
    def test_born_against_feynman_diagrams(self):
        from ggaa.born import born, helicity_amplitudes
        for b, c in ((0.1, 1 / 3), (0.45, -0.7), (0.8, 0.2)):
            h = helicity_amplitudes(b, c)
            self.assertAlmostEqual(abs(h[1, 1]) / abs(complex(born(1, b, c))), 1, places=10)
            self.assertAlmostEqual(abs(h[1, -1]) / abs(complex(born(2, b, c))), 1, places=10)

    def test_threshold_behaviour(self):
        from ggaa.born import born
        with mp.workdps(30):
            b, c = mp.mpf("1e-4"), mp.mpf(1) / 3
            self.assertLess(abs(born(1, b, c) - 2j * 125**2) / (2 * 125**2), 1e-7)        # i s/2 -> 2 i m^2
            ratio = born(2, b, c) / (2j * 125**2 * b**2 * (1 - c**2))                     # M_2 ~ 2 i m^2 beta^2 sin^2
            self.assertLess(abs(ratio - 1), 1e-3)


@unittest.skipUnless((INPUTS / "Finite_1Loop_P1").exists(), "run inputs/fetch_inputs.sh first")
class OneLoop(unittest.TestCase):
    def test_reference_values(self):
        from ggaa.oneloop import OneLoop as Evaluator
        amp = Evaluator()
        with mp.workdps(60):
            value = amp(1, mp.mpf(1) / 5, 0, scale_fix=False)
            expected = mp.mpc("-2531583.169999240279096260950887850230295",
                              "-2978915.948647069844502666628586329640816")
            self.assertLess(relative_error(value, expected), 1e-35)

    def test_scale_invariance(self):
        """With the scale term the remainder depends on m_A and mu_R only through their ratio."""
        from ggaa.oneloop import OneLoop as Evaluator
        from ggaa.born import born
        amp = Evaluator()
        with mp.workdps(50):
            b, c = mp.mpf(1) / 10, mp.mpf(1) / 3
            r = [amp(1, b, c, mA=m, muR=m) / m**2 for m in (125, 250)]
            self.assertLess(abs(r[1] - r[0]) / abs(r[0]), 1e-40)
            raw = [amp(1, b, c, mA=m, muR=m, scale_fix=False) / m**2 for m in (125, 250)]
            shift = -2 * mp.mpf(23) / 3 * mp.log(4) * born(1, b, c, 125) / 125**2
            self.assertLess(abs((raw[1] - raw[0]) - shift), 1e-40)


if __name__ == "__main__":
    unittest.main()


class HiggsPairs(unittest.TestCase):
    def test_forward_weight(self):
        """One-loop gg -> HH: the forward logarithm follows F_S (data/dps_hh, scripts/check_dps_hh.py)."""
        import csv
        path = pathlib.Path(__file__).resolve().parents[1] / "data" / "dps_hh" / "dps_hh.csv"
        for r in csv.DictReader(open(path)):
            x = (1 + float(r["beta"])) / 2
            fs = 1 / (x * (1 - x)) + 2 * x * (1 - x)
            fa = 1 / (x * (1 - x)) - 2 * x * (1 - x)
            self.assertLess(abs(-float(r["HH"]) / fs - 1), 1e-12)
            self.assertLess(abs(float(r["AA"]) / fa - 1), 1e-12)
            self.assertLess(float(r["HH_M2/HH_M1"]), 1e-10)


class Renormalisation2019(unittest.TestCase):
    def test_earlier_run_and_pole_record(self):
        """App. B: the earlier run's Z2 is off by 2 beta1/eps; the June 2019 record fixes FLAG = 1/96."""
        import subprocess
        script = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "check_renorm_2019.py"
        run = subprocess.run([sys.executable, str(script)], capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)

