"""Compute the best possible error_2 achievable by a quantized dictionary U_q by projecting G onto P(U_q)."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
import numpy as np
import torch
from compress_chem import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def check_dict(name, U, tree, brots=[8]):
    print(f"\n--- Checking {name} (shape {U.shape}) ---")
    sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
    ii, jj, fac = sp["ii"], sp["jj"], sp["fac"]
    ev, V = sp["eigenvalues"], sp["eigenvectors"]
    G = (V * np.maximum(ev, 0)) @ V.T
    total_norm_sq = float(np.square(np.maximum(ev, 0)).sum())
    
    # Also float projection floor:
    u_norm = U.reshape(-1, 76) / np.linalg.norm(U.reshape(-1, 76), axis=1, keepdims=True)
    P_float = u_norm[:, ii] * u_norm[:, jj] * fac
    # SVD or eigh of P P^T
    S = P_float @ P_float.T
    se, sv = np.linalg.eigh(S)
    pos = se > 1e-12
    sv = sv[:, pos]
    se = se[pos]
    iroot = (sv / np.sqrt(se)) @ sv.T
    Q = iroot @ P_float
    H = Q @ G @ Q.T
    inside_sq = float(np.square(H).sum())
    float_floor = np.sqrt(max(0.0, total_norm_sq - inside_sq))
    print(f"Float projection floor: {float_floor*1000:7.3f} mHa")
    
    for brot in brots:
        words = encode(U.reshape(-1, 76), tree, brot)
        U_q = reconstruct_u_vectors(words / (1 << (brot - 1)), tree)
        u_norm_q = U_q / np.linalg.norm(U_q, axis=1, keepdims=True)
        P_q = u_norm_q[:, ii] * u_norm_q[:, jj] * fac
        S = P_q @ P_q.T
        se, sv = np.linalg.eigh(S)
        pos = se > 1e-12
        sv = sv[:, pos]
        se = se[pos]
        iroot = (sv / np.sqrt(se)) @ sv.T
        Q = iroot @ P_q
        H = Q @ G @ Q.T
        inside_sq = float(np.square(H).sum())
        q_floor = np.sqrt(max(0.0, total_norm_sq - inside_sq))
        print(f"brot={brot:2d} quantized floor: {q_floor*1000:7.3f} mHa")

tree = np.load(ROOT / "initial_experiments/3/candidates/chem35c45_rot.npz")["tree_edges"]

# 1. global490
z1 = np.load(ROOT / "initial_experiments/3/candidates/global490_anneal_fit.npz")
check_dict("global490_anneal_fit", z1["U"], tree)

# 2. residual719 (the baseline dictionary!)
z2 = np.load(ROOT / "initial_experiments/4/agy/dfthcsos/submission.npz")
check_dict("submission719", z2["U"], tree)
