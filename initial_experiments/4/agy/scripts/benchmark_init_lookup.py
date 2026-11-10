"""Benchmark Kickmix's built-in init_lookup and del_lookup on 2048-entry table."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km

b = km.CircuitBuilder()
addr = b.create_quantum_register(11, 'addr')
target = b.create_quantum_register(28, 't')

table = list(range(2048))

b.init_lookup(table=table, address=addr, target=target)
b.del_lookup(table=table, address=addr, target=target)

c = b.finish_circuit()
print(f"init_lookup + del_lookup on 2048 entries:")
print(f"Total qubits: {c.num_qubits}")
print(f"Toffoli-like: {c.operation_statistics().num_toffoli_like}")
print(f"Operations: {c.operation_statistics().operation_counts}")
