"""Compute exact cost and ratio for various R, B, C with brot=8."""

import math
from dfthcsos.verify import heuristic_costs

N = 76
brot = 8
J = 13
K = math.ceil((N - 1) / J)
base = 3355599588.3155375
lam = 48.895

configs = [
    # Baseline
    (1, 719, 350),
    # R=2
    (2, 300, 200),
    (2, 250, 180),
    (2, 200, 175),
    (2, 150, 175),
    # R=3
    (3, 200, 120),
    (3, 150, 120),
    # R=4
    (4, 256, 128),
    (4, 150, 100),
    (4, 100, 90),
]

print(f"{'R':>2} {'B':>4} {'C':>4} | {'T_ang':>6} {'T_heur':>6} {'T_tot':>6} {'Q_tot':>5} | {'Cost':>12} {'Ratio':>7}")
print("-" * 65)

for R, B, C in configs:
    A = N + R * B
    sel = math.ceil(math.log2(A))
    # ANGLES QROM cost: 2 * K * (A - 2) Toffolis
    ang_toff = 2 * K * (A - 2)
    h = heuristic_costs(N, R, B, C, brot)
    T_tot = ang_toff + h['toffoli']
    Q_tot = h['persistent_qubits'] + (J * brot) + max(h['temporary_qubits'], sel - 1)
    cost = lam * T_tot * Q_tot
    ratio = base / cost
    print(f"{R:2d} {B:4d} {C:4d} | {ang_toff:6d} {h['toffoli']:6d} {T_tot:6d} {Q_tot:5d} | {cost:12,.0f} {ratio:6.2f}x")
