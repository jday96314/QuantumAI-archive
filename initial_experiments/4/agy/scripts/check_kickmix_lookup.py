"""Check b.del_lookup and b.lookup in kickmix."""
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km

b = km.CircuitBuilder()
addr = b.create_quantum_register(4, 'a')
target = b.create_quantum_register(8, 't')

# Check available lookup methods
methods = [m for m in dir(b) if 'lookup' in m or 'table' in m or 'qrom' in m]
print("Lookup methods in CircuitBuilder:", methods)
