"""Build linear invsqrt with defer_workspace and borrow_integration to minimize qubits."""

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

for qlim in [100, 90, 80, 75, 70, 65, 60, 55, 50]:
    try:
        c = polynomial.build(
            10, (28, 18), bias=0, rounded=True, booth=True,
            defer_workspace=True, borrow_integration=True, carry_bits=1,
            reuse_normalized_top=True, qubit_limit=qlim
        )
        print(f"qubit_limit={qlim:2d} -> SUCCESS! Allocated: {c.num_qubits}, Toffoli: {c.operation_statistics().num_toffoli_like}")
    except Exception as e:
        print(f"qubit_limit={qlim:2d} -> FAILED: {e}")
