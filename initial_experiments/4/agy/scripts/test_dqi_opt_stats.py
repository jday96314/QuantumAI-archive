"""Check qubits and Toffolis of build_opt_dqi_circuit in dqi_opt.py."""

import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers/dqi")
import kickmix as km
from dqi_opt import build_opt_dqi_circuit

print("Building opt dqi circuit...")
c = build_opt_dqi_circuit(n_iters=216)
print("Opt DQI circuit:")
print(f"Num qubits: {c.num_qubits}")
print(f"Num toffolis: {c.operation_statistics().num_toffoli_like:,}")
