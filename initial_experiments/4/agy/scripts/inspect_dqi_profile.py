import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km

c = km.Circuit.from_file("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/dqi.kmb", format="kmb")
print("Total qubits:", c.num_qubits)
print("Total bits:", c.num_bits)
print("Total Toffoli-like:", c.operation_statistics().num_toffoli_like)
print("Operation counts:", c.operation_statistics().operation_counts)
