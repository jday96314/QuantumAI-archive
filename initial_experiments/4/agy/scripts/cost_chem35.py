#!/usr/bin/env python3
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
from dfthcsos.verify import heuristic_costs
import math

base = 3355599588.3155375
N = 76
lam = 48.256
R, B, C = 35, 60, 45

for brot in [8, 10, 12]:
    h = heuristic_costs(N, R, B, C, brot)
    J = 10
    K = brot
    A = N + R * B
    sel = math.ceil(math.log2(A))
    ang_toff = 2 * K * (A - 2)
    T = ang_toff + h["toffoli"]
    Q = h["persistent_qubits"] + (J * brot) + max(h["temporary_qubits"], sel - 1)
    cost = lam * T * Q
    ratio = base / cost
    print(f"chem35 (brot={brot:2d}): T={T:,}, Q={Q:,} -> Cost: {cost:12,.0f} ({ratio:.3f}x)")
