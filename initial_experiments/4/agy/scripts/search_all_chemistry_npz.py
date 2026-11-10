#!/usr/bin/env python3
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import glob
import numpy as np
from dfthcsos.verify import heuristic_costs
import math

base = 3355599588.3155375
N = 76

files = glob.glob(str(ROOT / "initial_experiments/3/candidates/*.npz"))
print(f"Scanning {len(files)} files...")

best = []
for p in files:
    try:
        z = np.load(p)
        err2 = z.get("err2")
        if err2 is None:
            continue
        err2 = float(err2) * 1000 # in mHa
        if err2 > 20.05:
            continue
        
        U = z.get("U")
        w = z.get("w", z.get("W"))
        if U is None or w is None:
            continue
        
        R, B, _ = U.shape
        C = w.shape[-1]
        brot = int(z.get("brot", 8))
        lam = float(z.get("lambda_eff", 48.90))
        
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
        
        best.append((cost, ratio, Path(p).name, R, B, C, brot, err2, lam))
    except Exception:
        pass

best.sort(key=lambda x: x[0])
print(f"\nFound {len(best)} candidates with err2 <= 20.05 mHa:")
print(f"{'Cost':13s} | {'Ratio':7s} | {'Filename':40s} | R  B    C   brot | err2     lam")
print("-" * 90)
for cost, ratio, name, R, B, C, brot, err2, lam in best[:30]:
    print(f"{cost:11,.0f} | {ratio:6.3f}x | {name:40s} | {R:2d} {B:4d} {C:4d} {brot:4d} | {err2:6.2f}mHa {lam:5.2f}")
