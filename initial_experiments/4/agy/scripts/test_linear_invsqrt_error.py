"""Test exact maximum error of linear interpolation with F=11 on mantissa [1, 2)."""
import numpy as np

# Mantissa x in [1, 2)
# With lookup_frac F = 11, interval width h = 2^-11 = 1/2048
for F in [11, 10, 9]:
    h = 2.0**(-F)
    num_intervals = 1 << F
    max_rel_err = 0.0
    for i in range(num_intervals):
        x0 = 1.0 + i * h
        x1 = x0 + h
        mid = 0.5 * (x0 + x1)
        half_w = 0.5 * h
        t0 = mid - half_w / np.sqrt(2)
        t1 = mid + half_w / np.sqrt(2)
        y0 = 1.0 / np.sqrt(t0)
        y1 = 1.0 / np.sqrt(t1)
        slope = (y1 - y0) / (t1 - t0)
        intercept = y0 - slope * t0
        xs = np.linspace(x0, x1, 101)
        true_vals = 1.0 / np.sqrt(xs)
        approx_vals = intercept + slope * xs
        err = np.abs(true_vals - approx_vals)
        rel_err = np.max(err / true_vals)
        if rel_err > max_rel_err:
            max_rel_err = rel_err
    tol = 2.0**(-24)
    print(f"F={F:2d} (intervals={num_intervals:4d}) | max_rel_err={max_rel_err:.6e} | tol={tol:.6e} | ratio={max_rel_err/tol:.4f}")
