"""Inspect the ownership claims in circuit.py."""

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

print("Peak qubits:", b.peak_qubits)
print("Ownership claims at peak:")
if hasattr(b, "_claims"):
    for k, v in b._claims.items():
        print(f"  {k}: {len(v)} wires")
