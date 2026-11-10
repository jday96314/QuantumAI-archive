"""Run compact_wires on linear_invsqrt_f10."""

import subprocess
from pathlib import Path
import kickmix as km

raw_kmx = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/runs/linear_invsqrt_f10.raw.kmx")
compact_kmx = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/runs/linear_invsqrt_f10.compact.kmx")
compact_kmb = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/runs/linear_invsqrt_f10.compact.kmb")

print("Exporting to raw.kmx...")
c = km.Circuit.from_file("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/runs/linear_invsqrt_f10.kmb", format="kmb")
print(f"Initial qubits: {c.num_qubits}")
c.to_file(str(raw_kmx), format="kmx")

print("Running compact_wires...")
compact_wires_bin = "/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/compact_wires"
subprocess.run([compact_wires_bin, str(raw_kmx), str(compact_kmx), "22"], check=True)

print("Loading compacted circuit...")
c_compact = km.Circuit.from_file(str(compact_kmx))
print(f"Compacted qubits: {c_compact.num_qubits}")
c_compact.to_file(str(compact_kmb), format="kmb")
print(f"Saved {compact_kmb}")
