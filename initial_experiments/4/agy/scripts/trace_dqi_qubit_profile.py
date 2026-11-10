"""Analyze exact qubit allocation timeline of dqi.kmb."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km

c = km.Circuit.from_file("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/dqi.kmb", format="kmb")
print(f"DQI circuit: qubits={c.num_qubits}, bits={c.num_bits}")

# Let's inspect operation counts
stats = c.operation_statistics()
print("Toffoli count:", stats.num_toffoli_like)
print("Total ops:", c.num_operations)
