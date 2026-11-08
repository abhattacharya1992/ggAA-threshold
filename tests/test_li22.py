"""Regression tests for cancellation-safe Li22 series summation and caching."""

import json
import pathlib
import sys
import tempfile
import unittest
from unittest.mock import patch

import mpmath as mp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src" / "python"))

from ggaa.li22 import CACHE_VERSION, li22, li22_integral, li22_series
from ggaa.twoloop import Evaluator, symbols


class Li22RegressionTests(unittest.TestCase):
    def test_cancelled_inner_sum(self):
        # At i=3, y + y**2/4 vanishes for y=-4, but the tail is nonzero.
        for digits in (30, 60, 80):
            with self.subTest(digits=digits), mp.workdps(digits):
                x, y = mp.mpf(3) / 16, mp.mpf(-4)
                expected = li22_integral(x, y)
                self.assertLess(abs(li22(x, y) / expected - 1), mp.mpf(10)**(8 - digits))

    def test_defining_sum(self):
        with mp.workdps(80):
            x, y = mp.mpf(3) / 16, mp.mpf(-4)
            inner, total = mp.mpf(0), mp.mpf(0)
            xp, yp = mp.mpf(1), mp.mpf(1)
            for i in range(1, 1601):
                xp *= x
                total += xp * inner / i**2
                yp *= y
                inner += yp / i**2
            self.assertLess(abs(li22(x, y) / total - 1), mp.mpf('1e-72'))

    def test_zero_arguments(self):
        with mp.workdps(60):
            for x, y in ((0, 0), (0, -4), ('0.8', 0)):
                with self.subTest(x=x, y=y):
                    self.assertEqual(li22(x, y), 0)

    def test_other_convergent_arguments(self):
        with mp.workdps(60):
            pairs = [(mp.mpf('0.3'), mp.mpf('0.4')),
                     (mp.mpf('-0.3'), mp.mpf('1.2')),
                     (mp.mpc('0.2', '0.1'), mp.mpc('-0.4', '0.2'))]
            for x, y in pairs:
                with self.subTest(x=x, y=y):
                    expected = li22_integral(x, y)
                    self.assertLess(abs(li22(x, y) / expected - 1), mp.mpf('1e-52'))

    def test_physical_arguments_through_evaluator(self):
        # The published Li22 arguments at beta=1/7 and cos(theta)=1/2.
        with mp.workdps(80):
            kin = symbols(mp.mpf(1) / 7, mp.mpf(1) / 2)
            value = Evaluator(kin).evaluate('Li22[-x*(x-z)/(1+x),(1+x)/(x-z)]')
            expected = li22_integral(mp.mpf(3) / 16, mp.mpf(-4))
            self.assertLess(abs(value / expected - 1), mp.mpf('1e-72'))

    def test_old_persistent_value_is_not_reused(self):
        with tempfile.TemporaryDirectory() as tmp, mp.workdps(60):
            store = pathlib.Path(tmp) / 'li22_store.json'
            args = (mp.mpf(3) / 16, mp.mpf(-4))
            old_key = f'{mp.mp.dps}:' + ','.join(mp.nstr(a, mp.mp.dps - 8) for a in args)
            store.write_text(json.dumps({old_key: '-0.03515625'}))
            expected = li22_integral(*args)
            value = Evaluator({}, store).function('Li22', args)
            self.assertLess(abs(value / expected - 1), mp.mpf('1e-52'))
            self.assertIn(f'{CACHE_VERSION}:{old_key}', json.loads(store.read_text()))
            with patch('ggaa.twoloop.li22', side_effect=AssertionError('new cache was not reused')):
                cached = Evaluator({}, store).function('Li22', args)
            self.assertLess(abs(cached / expected - 1), mp.mpf('1e-52'))

    def test_nonconvergent_series_is_rejected(self):
        with self.assertRaises(ValueError):
            li22_series(mp.mpf(1), mp.mpf('0.5'))
        with self.assertRaises(ValueError):
            li22_series(mp.mpf('0.5'), mp.mpf(2))


if __name__ == '__main__':
    unittest.main()
