"""Inspect the exact peak live wires and allocation in dqi.kmb."""

import sys
from pathlib import Path
from collections import defaultdict
import json

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km

dqi_path = ROOT / "initial_experiments/4/agy/dqi.kmb"
print(f"Loading {dqi_path}...")
c = km.Circuit.from_file(dqi_path, format="kmb")
print(f"Total qubits: {c.num_qubits}, Total ops: {c.num_operations}")

# Let's inspect which qubits are touched most frequently
qubit_touches = [0] * c.num_qubits
# Scan circuit lines or operations if possible
stats = c.operation_statistics()
print(f"Toffolis: {stats.num_toffoli_like}")
