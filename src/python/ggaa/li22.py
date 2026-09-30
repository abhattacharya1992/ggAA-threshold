"""The weight-four multiple polylogarithm

    Li_{2,2}(x, y) = sum_{i > j >= 1} x^i y^j / (i^2 j^2),

in the convention of Frellesvig, Tommasini and Wever, JHEP 03 (2016) 189, which is the
one used by the ancillary files of arXiv:1909.08993.

For max(|x|, |x y|) < 0.99 the defining series is summed directly.  Closer to the edge
of convergence, which happens near threshold, the integral representation

    Li_{2,2}(x, y) = -x * Integral_0^1 dt ln(t) Li_2(x y t) / (1 - x t)

is used; it follows from the definition by writing 1/i^2 = -Integral_0^1 t^(i-1) ln t dt.
"""

import mpmath as mp

SERIES_LIMIT = 0.99


def li22_series(x, y):
    rate = max(abs(x), abs(x * y))
    if rate >= 1:
        raise ValueError(f"Li22 series does not converge at ({x}, {y})")
    tolerance = mp.eps * (1 - rate) / 4
    inner = mp.mpf(0)          # sum_{j < i} y^j / j^2
    x_power, y_power = mp.mpf(1), mp.mpf(1)
    total = mp.mpf(0)
    i = 1
    while True:
        x_power *= x
        term = x_power * inner / i**2
        total += term
        y_power *= y
        inner += y_power / i**2
        if i > 2 and abs(term) <= tolerance * abs(total):
            # the next terms are bounded by a geometric series in `rate`
            return total
        i += 1


def li22_integral(x, y):
    w = x * y
    integrand = lambda t: -mp.log(t) * mp.polylog(2, w * t) / (1 - x * t)
    return x * mp.quad(integrand, [0, 0.5, 0.9, 0.99, 1])


def li22(x, y):
    x, y = mp.mpmathify(x), mp.mpmathify(y)
    if max(abs(x), abs(x * y)) < SERIES_LIMIT:
        return li22_series(x, y)
    return li22_integral(x, y)
