import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
from dfthcsos.verify import heuristic_costs
import math

N = 76
brot = 8
J = 13
K = math.ceil((N - 1) / J)
base = 3355599588.3155375

configs = [
    (2, 128, 128),
    (2, 120, 150),
    (2, 100, 150),
    (2, 128, 100),
    (1, 250, 200),
    (1, 250, 150),
    (1, 200, 200),
    (1, 200, 150),
    (2, 150, 150),
]

print(f"{'R':>2} {'B':>4} {'C':>4} {'RxC':>5} {'RxB':>5} {'T_ang':>6} {'T_heur':>6} {'T_tot':>6} {'Q_tot':>5} {'Cost':>12} {'Ratio':>7}")
print("-" * 75)

for R, B, C in configs:
    A = N + R * B
    sel = math.ceil(math.log2(A))
    ang_toff = 2 * K * (A - 2)
    try:
        h = heuristic_costs(N, R, B, C, brot)
        T_heur = h['toffoli']
        T = ang_toff + T_heur
        Q = h['persistent_qubits'] + (J * brot) + max(h['temporary_qubits'], sel - 1)
        cost = 48.895 * T * Q
        ratio = base / cost
        print(f"{R:2d} {B:4d} {C:4d} {R*C:5d} {R*B:5d} {ang_toff:6d} {T_heur:6d} {T:6d} {Q:5d} {cost:12,.0f} {ratio:6.2f}x")
    except Exception as e:
        print(f"R={R}, B={B}, C={C}: {e}")
