"""Scan all candidate files and compute exact host verifier metrics and score."""

import glob, os, sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
import numpy as np
from dfthcsos.verify import load_inputs, evaluate_words, heuristic_costs

ref_dict = np.load(ROOT / "initial_experiments/3/parameters/chem_reference.npz")
base = dict(np.load(ROOT / "initial_experiments/3/candidates/chem35c45_rot.npz"))

files = sorted(glob.glob("../../3/candidates/*.npz"))
print(f"Scanning {len(files)} files for valid passing candidates...")

incumbent_score = 757163489.50
passing = []

for f in files:
    try:
        z = np.load(f)
        if "onebody_angle_words" not in z or "angle_words" not in z:
            continue
        words1 = z["onebody_angle_words"]
        words2 = z["angle_words"].reshape(-1, 75)
        words = np.r_[words1, words2]
        tree = z.get("tree_edges", base["tree_edges"])
        brot = int(z.get("brot", 8))
        
        U = z["U"]
        w = z["w"]
        W = z["W"]
        R, B, N = U.shape if U.ndim == 3 else (1, U.shape[0], U.shape[1])
        C = w.shape[-1]
        
        # Fast pre-check: err2
        err2_pre = float(z.get("err2", 1.0)) * 1000
        if err2_pre > 20.05:
            continue
            
        fact = dict(U=U, W=w, H2_I=W + w.sum(1), H1_I=base["H1_I"], Hsym_lower=base["Hsym_lower"])
        inputs = load_inputs(fact)
        metrics = evaluate_words(inputs, words, tree, brot)
        
        err1 = metrics["error_1"] * 1000
        err2 = metrics["error_2"] * 1000
        gap = metrics["E_gap"]
        lam_eff = metrics["lambda_eff"]
        in_domain = metrics["gap_in_domain"]
        
        if err1 <= 20.0 and err2 <= 20.0 and in_domain and np.isfinite(lam_eff):
            # Compute score
            h = heuristic_costs(N, R, B, C, brot)
            A = N + R * B
            sel = int(np.ceil(np.log2(A)))
            ang_toff = 2 * brot * (A - 2) # approx
            T = ang_toff + h["toffoli"]
            Q = h["persistent_qubits"] + (13 * brot) + max(h["temporary_qubits"], sel - 1)
            score = lam_eff * T * Q
            ratio = 3355599588.3155375 / score
            passing.append((score, ratio, B, C, R, brot, err1, err2, lam_eff, os.path.basename(f)))
            print(f"*** PASS: {os.path.basename(f)} | Score = {score:,.0f} | Ratio = {ratio:.4f}x | B={B}, C={C}, lam={lam_eff:.2f}, err2={err2:.2f} mHa ***")
    except Exception as e:
        pass

passing.sort(key=lambda x: x[0])
print(f"\nTotal passing candidates found: {len(passing)}")
for p in passing[:10]:
    print(p)
