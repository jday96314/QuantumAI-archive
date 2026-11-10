import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "initial_experiments/3"),
    str(ROOT / "initial_experiments/1"),
]
import kickmix as km
import binary_centered_invsqrt as polynomial
import numpy as np

table = np.load(ROOT / "initial_experiments/3/candidates/margin_curve8_t5.9e-08.npy")
polynomial.table_for = lambda *args: table

for qlim in [94, 92, 90, 88, 86, 84, 82, 80, 78]:
    try:
        c = polynomial.build(
            7, (27, 18, 8), lookup_frac=8, accum_bits=[18, 9, None],
            carry_bits=1, rounded=True, booth=True, reuse_normalized_top=True,
            output_extra_bits=10, defer_workspace=True, qubit_limit=qlim,
            curve_keep_bits=4, defer_curve_phase=True, product_limit=17,
            borrow_integration=True
        )
        print(f"qubit_limit={qlim:2d} -> SUCCEEDED! Circuit num_qubits: {c.num_qubits}, Toffolis: {c.operation_statistics().num_toffoli_like}")
    except Exception as e:
        print(f"qubit_limit={qlim:2d} -> FAILED: {e}")
