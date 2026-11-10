"""Test pruning the verified residual719 factorization from B=719 down to B=350-500."""

import sys
from pathlib import Path
import numpy as np

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
from dfthcsos.scoreinfra.reference import load_reference, PACKAGED_FCIDUMP
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
from dfthcsos.verify import load_inputs, evaluate_words, heuristic_costs

z = dict(np.load(ROOT / "initial_experiments/4/agy/dfthcsos/submission.npz"))
U = z["U"][0] # (719, 76)
w = z["W"][0] # (719, 350)
H2_I = z["H2_I"][0] # (350,)
tree = z["tree_edges"]
words = z["angle_words"] # (795, 75)
brot = 8

print(f"Original residual719: B={len(U)}, C={w.shape[1]}")

# Compute importance of each basis vector b:
# The contribution to X_pq is sum_b w_bc u_bp u_bq
# So vector b has weight ||w[b, :]||_2
weights = np.linalg.norm(w, axis=1)
order = np.argsort(weights)[::-1]

base = 3355599588.3155375
N = 76

# Test keeping top B_keep vectors
for B_keep in [650, 600, 550, 500, 450, 400, 350]:
    keep_idx = sorted(order[:B_keep])
    U_sub = U[keep_idx]
    w_sub = w[keep_idx]
    
    # One-body words: first 76
    one_words = words[:76]
    # Two-body words: 76 + keep_idx
    two_words = words[76 + np.array(keep_idx)]
    sub_words = np.r_[one_words, two_words]
    
    fact = dict(
        U=U_sub[None, :, :],
        W=w_sub[None, :, :],
        H2_I=H2_I[None, :],
        H1_I=z["H1_I"],
        Hsym_lower=z["Hsym_lower"],
    )
    
    try:
        inputs = load_inputs(fact)
        metrics = evaluate_words(inputs, sub_words, tree, brot)
        err1 = metrics["error_1"] * 1000
        err2 = metrics["error_2"] * 1000
        gap = metrics["E_gap"]
        lam_eff = metrics["lambda_eff"]
        
        # Calculate cost
        R = 1
        C = w_sub.shape[1]
        J = 13
        K = 6
        A = N + R * B_keep
        sel = (A - 1).bit_length()
        ang_toff = 2 * K * (A - 2)
        h = heuristic_costs(N, R, B_keep, C, brot)
        T = ang_toff + h["toffoli"]
        Q = h["persistent_qubits"] + (J * brot) + max(h["temporary_qubits"], sel - 1)
        cost = lam_eff * T * Q
        ratio = base / cost
        
        passed = (err1 < 20.0 and err2 < 20.0 and gap > 0 and metrics["gap_in_domain"])
        print(f"B={B_keep:3d} -> err1={err1:5.2f} mHa, err2={err2:5.2f} mHa, gap={gap:5.2f} Ha, lam={lam_eff:5.2f} | Cost: {cost:11,.0f} ({ratio:.2f}x) | PASS: {passed}")
    except Exception as e:
        print(f"B={B_keep:3d} -> Failed: {e}")
