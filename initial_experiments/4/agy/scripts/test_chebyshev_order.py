"""Test Chebyshev interpolation error for various F and order."""

import numpy as np

tol = 2.0**(-24)

for order in [1, 2]:
    for frac in [6, 7, 8, 9, 10]:
        h = 2.0**(-frac)
        # Check worst case theoretical Chebyshev error over domain [1, 2]
        # In interval [a - h/2, a + h/2], 1/sqrt(x)
        # For order p Chebyshev, max error <= (h/2)^(p+1) / (2^p * (p+1)!) * max |f^(p+1)|
        # Let's compute actual Chebyshev interpolation for all intervals
        max_err = 0.0
        intervals = 1 << frac
        for i in range(intervals):
            x_mid = 1.0 + (i + 0.5) * h
            domain = [x_mid - h/2, x_mid + h/2]
            c = np.polynomial.Chebyshev.interpolate(lambda x: 1.0 / np.sqrt(x), order, domain=domain)
            # Sample densely in this interval
            xs = np.linspace(domain[0], domain[1], 100)
            err = np.max(np.abs(c(xs) - 1.0 / np.sqrt(xs)))
            rel_err = err / (1.0 / np.sqrt(xs[-1])) # relative error
            if rel_err > max_err:
                max_err = rel_err
        
        ratio = max_err / tol
        print(f"order={order}, frac={frac:2d} -> max relative error = {max_err:.3e} ({ratio:7.3f}x tolerance)")
