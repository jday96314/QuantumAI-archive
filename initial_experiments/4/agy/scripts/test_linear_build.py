"""Test binary_centered_invsqrt with linear widths (len(widths) == 2)."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
    str(ROOT / "initial_experiments/1"),
]
import kickmix as km
import binary_centered_invsqrt as polynomial
import numpy as np

# Test model error with frac=10, widths=(27, 18)
xs = np.arange(1, 3133453)
kw = dict(frac=10, widths=(27, 18), rounded=True, booth=True)
print("Evaluating linear model on all 3,133,452 inputs...")
res = polynomial.model(xs, **kw)

# Compute error relative to 1/sqrt(x)
error = (res / (2.0**(27 + 10))) * np.sqrt(xs) - 1.0
max_rel = float(np.max(np.abs(error)))
tol = 2.0**(-24)
print(f"Max relative error: {max_rel:.6e}")
print(f"Required tolerance: {tol:.6e}")
print(f"Ratio max_rel / tol: {max_rel / tol:.4f}")

if max_rel < tol:
    print("SUCCESS! Model passes all 3.13M points!")
    print("Now testing circuit build...")
    c = polynomial.build(10, (27, 18), rounded=True, booth=True)
    print(f"Circuit built! Qubits: {c.num_qubits}, Toffolis: {c.operation_statistics().num_toffoli_like}")
else:
    print(f"Model failed tolerance check (ratio {max_rel / tol:.2f})")
