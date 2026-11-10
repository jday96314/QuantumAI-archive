#!/usr/bin/env python3
"""Test encoding and evaluating positive650_continue.npz."""

import sys
from pathlib import Path
import numpy as np

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
from compress_chem import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
from dfthcsos.verify import load_inputs, evaluate_words, heuristic_costs

fit_path = ROOT / "initial_experiments/3/candidates/positive650_continue.npz"
z = dict(np.load(fit_path))
U = z["U"]
w = z["w"]
W = z["W"]
R, B, N = U.shape
C = w.shape[2]
print(f"Loaded positive650: R={R}, B={B}, C={C}, float err2={float(z['err2'])*1000:.3f} mHa")

base = dict(np.load(ROOT / "initial_experiments/3/candidates/chem35c45_rot.npz"))
tree = z.get("tree_edges", base["tree_edges"])

fact = dict(U=U, W=w, H2_I=W + w.sum(1), H1_I=base["H1_I"], Hsym_lower=base["Hsym_lower"])
inputs = load_inputs(fact)

ref_dict = np.load(ROOT / "initial_experiments/3/parameters/chem_reference.npz")
h0 = ref_dict['h1_thc'] - 0.5 * int(ref_dict['electrons']) * ref_dict['S'] + float(ref_dict['beta_1']) * np.eye(N)

for brot in [8, 9, 10, 11, 12]:
    two = encode(U.reshape(-1, N), tree, brot)
    u_rec = reconstruct_u_vectors(two / (1 << (brot - 1)), tree).reshape(R, B, N)
    X = np.einsum('rbp,rbq,rbc->rcpq', u_rec, u_rec, w, optimize=True)
    h = h0 - np.einsum('rc,rcpq->pq', W, X, optimize=True)
    _, v = np.linalg.eigh(0.5 * (h + h.T))
    one = encode(v.T, tree, brot)
    
    words = np.r_[one, two.reshape(R * B, 75)]
    
    try:
        metrics = evaluate_words(inputs, words, tree, brot)
        err1 = metrics["error_1"] * 1000
        err2 = metrics["error_2"] * 1000
        gap = metrics["E_gap"]
        lam_eff = metrics["lambda_eff"]
        passed = (err1 < 20.0 and err2 < 20.0 and gap > 0 and metrics["gap_in_domain"])
        print(f"brot={brot:2d} -> err1: {err1:6.2f} mHa, err2: {err2:6.2f} mHa, gap: {gap:6.3f} Ha, lam_eff: {lam_eff:6.2f} | PASS: {passed}")
    except Exception as e:
        print(f"brot={brot:2d} -> Failed: {e}")
