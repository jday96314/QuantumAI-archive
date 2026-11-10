"""Sweep chemistry cost formula across R, B, C, brot to find the theoretical minimum costs."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import math
from dfthcsos.verify import heuristic_costs

N = 76
base = 3355599588.3155375
lam = 48.895

print(f"{'R':>2} {'B':>3} {'C':>3} {'brot':>4} | {'T_ang':>6} {'T_prep':>6} {'T_giv':>5} {'T_tot':>6} | {'Q':>4} | {'Cost':>12} | {'Ratio':>6}")
print("-" * 75)

candidates = [
    (1, 719, 350, 8),
    (1, 719, 320, 8),
    (1, 719, 280, 8),
    (1, 719, 250, 8),
    (1, 719, 200, 8),
    (1, 500, 500, 8),
    (1, 490, 490, 8),
    (1, 400, 400, 8),
    (1, 350, 350, 8),
    (1, 300, 300, 8),
    (4, 256, 128, 8),
    (4, 256, 64, 8),
    (2, 256, 128, 8),
    (1, 256, 128, 8),
]

for R, B, C, brot in candidates:
    J = 13 if brot <= 8 else 10
    K = math.ceil((N - 1) / J)
    A = N + R * B
    sel = math.ceil(math.log2(A))
    ang_toff = 2 * K * (A - 2)
    h = heuristic_costs(N, R, B, C, brot)
    T = ang_toff + h["toffoli"]
    Q = h["persistent_qubits"] + (J * brot) + max(h["temporary_qubits"], sel - 1)
    cost = lam * T * Q
    ratio = base / cost
    print(f"{R:2d} {B:3d} {C:3d} {brot:4d} | {ang_toff:6d} {h['toffoli']:6d} {2*(N-1)*(2*brot-1):5d} {T:6d} | {Q:4d} | {cost:12,.0f} | {ratio:6.2f}x")
