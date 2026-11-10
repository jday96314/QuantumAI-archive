import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3")
sys.path[:0] = [str(ROOT), str(ROOT / "vendor/uncompute")]

import import_published as p
m = p.load()

# Let's inspect BM object and allocator
b = m.Builder()
syn = list(m.declare_registers(b, [(f's{i:02d}', 12) for i in range(70)]).values())

class Alloc(m.Alloc):
    def __init__(self, b, syn):
        super().__init__(b)
        self.syn = syn
        self.bb_count = 0
        self.borrowed = set()
        self.alloc_log = {}

    def open(self, name, n):
        qs = super().open(name, n)
        self.alloc_log[name] = self.alloc_log.get(name, 0) + len(qs)
        return qs

A = Alloc(b, syn)
held = b.qubits(23)
Dw = held[:12]
dW = held[12:20]
E = m.Em(b, *held[20:23], A)

# Print initial live qubits
print("Initially allocated:", b.next_qubit)
