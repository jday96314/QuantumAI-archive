#!/usr/bin/env python3
"""Profile exact Toffoli breakdown of _forney_kernel in stream_forney.py."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "initial_experiments/3"),
    str(ROOT / "initial_experiments/2"),
]
import kickmix as km
from fast_dqi import FastBuilder
from stream_forney import K, GF

b = FastBuilder(km.CircuitBuilder())
acc_s = b.create_quantum_register(12, "acc_s")
acc_d = b.create_quantum_register(12, "acc_d")
acc_o = b.create_quantum_register(12, "acc_o")
mask = b.create_classical_register(12, "mask")

K._forney_kernel(b, acc_s, acc_d, acc_o, mask)
c = b.finish_circuit()

stats = c.operation_statistics()
print("Total Toffoli per point:", stats.num_toffoli_like)
for k, v in stats.items():
    print(f"  {k}: {v}")
