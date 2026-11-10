"""Check stats of adaptivebooth_f7_m8_w27_18_9.kmb"""

from pathlib import Path
import sys
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km

c_path = ROOT / "initial_experiments/3/candidates/adaptivebooth_f7_m8_w27_18_9.kmb"
c = km.Circuit.from_file(c_path, format="kmb")
stats = c.operation_statistics()
print("Qubits:", c.num_qubits)
print("Bits:", c.num_bits)
print("Num ops:", stats.num_operations)
print("Toffoli like:", stats.num_toffoli_like)
print("Rotation count:", stats.rotation_count)
