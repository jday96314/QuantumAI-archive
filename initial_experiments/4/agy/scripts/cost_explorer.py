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
R = 1
brot = 8
J = 13
K = math.ceil((N - 1) / J)

base = 3355599588.3155375

print(f"{'B':>4} {'C':>4} {'T_ang':>6} {'T_heur':>6} {'T_tot':>6} {'Q_tot':>5} {'Cost':>12} {'Ratio':>7}")
print("-" * 55)

for B in [719, 500, 375, 350, 250, 200, 150, 100, 76]:
    for C in [350, 300, 250, 200, 150, 100]:
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
            print(f"{B:4d} {C:4d} {ang_toff:6d} {T_heur:6d} {T:6d} {Q:5d} {cost:12,.0f} {ratio:6.2f}x")
        except Exception as e:
            # print(f"B={B}, C={C}: {e}")
            pass
