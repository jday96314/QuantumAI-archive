"""Score all chemistry candidates in initial_experiments/3/candidates and report any that beat 757M."""

import sys
import glob
from pathlib import Path
import numpy as np

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
from dfthcsos.verify import heuristic_costs
import math

base = 3355599588.3155375
best_cost = 757163489.50

files = glob.glob(str(ROOT / "initial_experiments/3/candidates/*.npz"))
print(f"Scanning {len(files)} npz files...")

scored = []
for p in files:
    try:
        z = np.load(p)
        if "U" not in z or "w" not in z:
            continue
        err2 = float(z["err2"]) if "err2" in z else None
        if err2 is None or err2 > 0.02005:
            continue
        
        U = z["U"]
        w = z["w"]
        R, B, N = U.shape
        _, _, C = w.shape
        brot = int(z.get("brot", 8))
        lam_eff = float(z["lambda_eff"]) if "lambda_eff" in z else 48.895
        
        J = 13 if brot <= 8 else 10
        K = math.ceil((N - 1) / J)
        A = N + R * B
        sel = math.ceil(math.log2(A))
        ang_toff = 2 * K * (A - 2)
        h = heuristic_costs(N, R, B, C, brot)
        T = ang_toff + h["toffoli"]
        Q = h["persistent_qubits"] + (J * brot) + max(h["temporary_qubits"], sel - 1)
        cost = lam_eff * T * Q
        ratio = base / cost
        scored.append((cost, ratio, p, R, B, C, brot, lam_eff, err2))
    except Exception:
        pass

scored.sort(key=lambda x: x[0])
print(f"\nTop 15 lowest-cost valid candidates:")
for cost, ratio, p, R, B, C, brot, lam_eff, err2 in scored[:15]:
    fname = Path(p).name
    print(f"Cost: {cost:12,.0f} | Ratio: {ratio:6.2f}x | R={R:2d}, B={B:3d}, C={C:3d}, brot={brot:2d} | lam={lam_eff:.2f}, err2={err2*1000:.2f}mHa | {fname}")
