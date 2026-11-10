"""Test compact_wires on dqi_published.kmb to see its true peak live qubits."""

import subprocess
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km

pub_kmb = ROOT / "initial_experiments/3/candidates/dqi_published.kmb"
pub_kmx = Path("/tmp/pub_temp.kmx")
compact_kmx = Path("/tmp/pub_compact.kmx")

print("Exporting dqi_published.kmb to kmx...")
c = km.Circuit.from_file(pub_kmb, format="kmb")
print(f"Original dqi_published: qubits={c.num_qubits}, bits={c.num_bits}")

c.to_file(str(pub_kmx), format="kmx")
print("Running compact_wires with pinned=840...")

compact_wires_bin = ROOT / "initial_experiments/3/compact_wires"
res = subprocess.run([str(compact_wires_bin), str(pub_kmx), str(compact_kmx), "840"], capture_output=True, text=True)
print("compact_wires stdout:", res.stdout)
print("compact_wires stderr:", res.stderr)
