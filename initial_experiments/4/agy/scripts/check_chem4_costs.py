"""Check cost of chem4b256c128 at brot=8, 10, 11, 12."""

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
R, B, C = 4, 256, 128

for brot in [8, 9, 10, 11, 12]:
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
    print(f"brot={brot:2d}: T_ang={ang_toff}, T_prep={h['toffoli']}, T_tot={T}, Q={Q} -> Cost: {cost:12,.0f} ({ratio:.2f}x)")
