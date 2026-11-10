"""Test running compact_wires on dqi.kmb to see how many qubits are actually needed at peak."""

import subprocess
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km

dqi_kmb = ROOT / "initial_experiments/4/agy/dqi.kmb"
dqi_kmx = Path("/tmp/dqi_temp.kmx")
compact_kmx = Path("/tmp/dqi_compact.kmx")

print("Exporting dqi.kmb to kmx...")
c = km.Circuit.from_file(dqi_kmb, format="kmb")
print(f"Original dqi: qubits={c.num_qubits}, bits={c.num_bits}")

c.to_file(str(dqi_kmx), format="kmx")
print("Running compact_wires with pinned=840 (70 syndromes)...")

compact_wires_bin = ROOT / "initial_experiments/3/compact_wires"
res = subprocess.run([str(compact_wires_bin), str(dqi_kmx), str(compact_kmx), "840"], capture_output=True, text=True)
print("compact_wires stdout:", res.stdout)
print("compact_wires stderr:", res.stderr)
