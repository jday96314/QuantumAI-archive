"""Check runs/linear_invsqrt_f10.kmb statistics."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km

c = km.Circuit.from_file(ROOT / "initial_experiments/4/agy/runs/linear_invsqrt_f10.kmb", format="kmb")
print("Qubits:", c.num_qubits)
print("Toffoli:", c.operation_statistics().num_toffoli_like)
print("Ops:", len(c.operations))
