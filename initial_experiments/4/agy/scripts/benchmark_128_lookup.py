"""Benchmark init_lookup + del_lookup for 128 entries."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km

b = km.CircuitBuilder()
addr = b.create_quantum_register(7, 'addr')
target = b.create_quantum_register(54, 't')

table = list(range(128))

b.init_lookup(table=table, address=addr, target=target)
b.del_lookup(table=table, address=addr, target=target)

c = b.finish_circuit()
print(f"128 entries (54 bits target):")
print(f"Total qubits: {c.num_qubits}")
print(f"Toffoli-like: {c.operation_statistics().num_toffoli_like}")
