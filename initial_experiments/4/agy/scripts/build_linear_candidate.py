"""Build and test the pure linear invsqrt circuit with frac=10, widths=(28, 18)."""

import sys
import time
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

print("Building circuit: frac=10, widths=(28, 18), bias=0, rounded=True, booth=True...")
start = time.time()
c = polynomial.build(10, (28, 18), bias=0, rounded=True, booth=True)
elapsed = time.time() - start

print(f"Circuit built in {elapsed:.2f}s!")
print(f"Total qubits: {c.num_qubits}")
print(f"Toffoli-like: {c.operation_statistics().num_toffoli_like}")

# Save circuit
out_path = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/runs/linear_invsqrt_f10.kmb")
out_path.parent.mkdir(parents=True, exist_ok=True)
c.to_file(str(out_path), format="kmb")
print(f"Saved to {out_path}")
