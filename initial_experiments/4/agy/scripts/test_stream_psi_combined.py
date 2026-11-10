"""Test combining Bank syndrome reuse from stream_capture.py with the full Psi route of circuit.py."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "initial_experiments/3/vendor/uncompute"),
]
import importlib.util

spec = importlib.util.spec_from_file_location(
    "discovered_circuit",
    ROOT / "initial_experiments/3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/circuit.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

OldBank = mod.BankReg
bank_offset = 3

class SyndBank(OldBank):
    __slots__ = ('slot',)
    def __init__(self, A, name, width=mod.B):
        super().__init__(A, name, width)
        self.slot = None
        if name == 'bb':
            self.slot = A.bb_count
            A.bb_count += 1
    def mat(self):
        if self.wires is None and self.slot is not None and self.slot >= bank_offset:
            A = self.A
            self.wires = list(A.syn[self.slot - bank_offset])
            A.borrowed.update(self.wires)
            self.round = A.round
            A.bank_log.setdefault(A.round, []).append(self)
            return self.wires
        return super().mat()

mod.BankReg = SyndBank

OldAlloc = mod.Alloc
class SyndAlloc(OldAlloc):
    def __init__(self, b):
        super().__init__(b)
        self.bb_count = 0
        self.borrowed = set()
        self.syn = None

mod.Alloc = SyndAlloc

OldBM = mod.BM
class SyndBM(OldBM):
    def __init__(self, E, A, syn, Dw, dW):
        A.syn = syn
        super().__init__(E, A, syn, Dw, dW)

mod.BM = SyndBM

print("Testing modified circuit.py build with SyndBank...")
try:
    mod.ZMASKI = [0] * 70
    b, t = mod.build()
    print("SUCCESS!")
    print(f"Total ops: {len(b._items)}")
    print(f"Total Toffoli: {t}")
    print(f"Peak qubits: {b.peak_qubits}")
    print(f"Next qubit: {b.next_qubit}")
except Exception as e:
    print("FAILED with error:", type(e), e)
    import traceback
    traceback.print_exc()
