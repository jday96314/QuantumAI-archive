#!/usr/bin/env python3
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import numpy as np
import glob
from dfthcsos.verify import load_inputs, evaluate_words, heuristic_costs

cands = glob.glob(str(ROOT / "initial_experiments/3/candidates/*_q8_*.npz"))
print(f"Testing {len(cands)} candidates...")

valid = []
for p in cands:
    try:
        z = dict(np.load(p))
        if "angle_words" not in z or "tree_edges" not in z:
            continue
        U = z["U"]
        W = z.get("W", z.get("w"))
        if W is None:
            continue
        R, B, N = U.shape
        C = W.shape[-1]
        tree = z["tree_edges"]
        words = z["angle_words"]
        
        # Load inputs
        base = dict(np.load(ROOT / "initial_experiments/3/candidates/chem35c45_rot.npz"))
        fact = dict(U=U, W=W, H2_I=z.get("H2_I", W), H1_I=z.get("H1_I", base["H1_I"]), Hsym_lower=z.get("Hsym_lower", base["Hsym_lower"]))
        inputs = load_inputs(fact)
        
        m = evaluate_words(inputs, words, tree, brot=8)
        err1 = m["error_1"] * 1000
        err2 = m["error_2"] * 1000
        gap = m["E_gap"]
        lam_eff = m["lambda_eff"]
        gap_ok = m["gap_in_domain"]
        
        if err1 <= 20.0 and err2 <= 20.0 and gap_ok:
            h = heuristic_costs(N, R, B, C, 8)
            J = 10
            K = 8
            A = N + R * B
            sel = (A - 1).bit_length()
            ang_toff = 2 * K * (A - 2)
            T = ang_toff + h["toffoli"]
            Q = h["persistent_qubits"] + (J * 8) + max(h["temporary_qubits"], sel - 1)
            cost = lam_eff * T * Q
            ratio = 3355599588.3155375 / cost
            valid.append((cost, ratio, p, R, B, C, err1, err2, gap, lam_eff))
            print(f"PASS: {Path(p).name} | Cost: {cost:11,.0f} ({ratio:6.3f}x) | R={R} B={B} C={C} | err1={err1:5.2f} err2={err2:5.2f} mHa")
    except Exception as e:
        pass

valid.sort(key=lambda x: x[0])
print(f"\nTotal valid: {len(valid)}")
for cost, ratio, p, R, B, C, err1, err2, gap, lam_eff in valid:
    print(f"{Path(p).name:40s} | Cost: {cost:11,.0f} ({ratio:6.3f}x) | err2: {err2:5.2f} mHa")
