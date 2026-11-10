"""Trace the exact lifetime and max co-live qubits of circuit.py stage by stage."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "initial_experiments/3/vendor/uncompute"),
]
import importlib.util

spec = importlib.util.spec_from_file_location(
    "discovered_circuit",
    ROOT / "initial_experiments/3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/circuit.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

mod.ZMASKI = [0] * 70
b, t = mod.build()

print(f"Build complete!")
print(f"Total ops: {len(b._items)}")
print(f"Total Toffoli: {t}")
print(f"Peak qubits: {b.peak_qubits}")
print(f"Next qubit: {b.next_qubit}")
print(f"Free qubits at end: {len(b.free_qubits)}")

# Let's inspect phase by phase
if hasattr(b, "_phase_marks"):
    print("Phase marks:", b._phase_marks)
