"""Find lowest cost configurations where Total Rank R * C >= 290."""

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

results = []

for R in [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 30, 35]:
    for B in range(16, 400, 8):
        for C in range(8, 360, 4):
            if R * C < 290:
                continue
            if B * R > 800:
                continue
            for brot in [8]:
                try:
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
                    results.append((cost, ratio, R, B, C, brot, T, Q))
                except Exception:
                    pass

results.sort(key=lambda x: x[0])
print(f"Top 25 configurations with R * C >= 290 (brot=8):")
print(f"{'R':>2} {'B':>3} {'C':>3} {'R*C':>4} | {'T_tot':>6} | {'Q':>4} | {'Cost':>12} | {'Ratio':>6}")
print("-" * 60)
for cost, ratio, R, B, C, brot, T, Q in results[:25]:
    print(f"{R:2d} {B:3d} {C:3d} {R*C:4d} | {T:6d} | {Q:4d} | {cost:12,.0f} | {ratio:6.2f}x")
