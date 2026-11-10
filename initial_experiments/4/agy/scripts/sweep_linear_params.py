"""Test parameter combinations for pure linear invsqrt."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
    str(ROOT / "initial_experiments/1"),
]
import binary_centered_invsqrt as polynomial
import numpy as np

xs = np.arange(1, 3133453)
tol = 2.0**(-24)

for frac in [10, 11]:
    for widths in [(27, 18), (27, 19), (28, 18), (28, 19), (27, 20), (28, 20)]:
        for bias in [0, 1, 2, 3]:
            kw = dict(frac=frac, widths=widths, bias=bias, rounded=True, booth=True)
            res = polynomial.model(xs, **kw)
            error = (res / (2.0**(widths[0] + 10))) * np.sqrt(xs) - 1.0
            max_rel = float(np.max(np.abs(error)))
            ratio = max_rel / tol
            if ratio < 1.0:
                print(f"SUCCESS: frac={frac}, widths={widths}, bias={bias} | ratio={ratio:.4f} (< 1.0!)")
                break
            elif ratio < 1.15:
                print(f"Close: frac={frac}, widths={widths}, bias={bias} | ratio={ratio:.4f}")
