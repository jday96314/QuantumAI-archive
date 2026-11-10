#!/usr/bin/env python3
"""Check exact subspace projection floor of subsets of residual719 vectors."""

import sys
from pathlib import Path
import numpy as np

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]

z = dict(np.load(ROOT / "initial_experiments/4/agy/dfthcsos/submission.npz"))
U = z["U"][0]  # (719, 76)
w = z["W"][0]  # (719, 350)
N = 76
B_total = len(U)

# Reconstruct exact normalized u vectors
U_norm = U / np.linalg.norm(U, axis=1, keepdims=True)

# Load target spectrum
sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
ii, jj, fac = sp["ii"], sp["jj"], sp["fac"]
ev, V = sp["eigenvalues"], sp["eigenvectors"]
G = (V * np.maximum(ev, 0)) @ V.T
total_norm_sq = float(np.square(np.maximum(ev, 0)).sum())

# Feature matrix for all 719 vectors: shape (719, 2926)
P_all = U_norm[:, ii] * U_norm[:, jj] * fac

# Compute vector importance:
# Method 1: norm of w row
w_norms = np.linalg.norm(w, axis=1)
order_w = np.argsort(w_norms)[::-1]

# Method 2: greedy leverage score / correlation with residual
print("Evaluating projection floor vs subset size B...")
print(f"Total 719 vectors: B={B_total}")

for B_keep in [719, 650, 600, 550, 500, 450, 400, 350]:
    idx = order_w[:B_keep]
    P_sub = P_all[idx] # (B_keep, 2926)
    
    # Projection floor of span(P_sub):
    # Orthogonal projection of G onto row span of P_sub
    # Gram matrix S = P_sub @ P_sub.T (B_keep x B_keep)
    S = P_sub @ P_sub.T
    se, sv = np.linalg.eigh(S)
    pos = se > 1e-10
    sv = sv[:, pos]
    se = se[pos]
    iroot = (sv / np.sqrt(se)) @ sv.T
    Q = iroot @ P_sub # (rank, 2926) orthonormal basis
    
    # Project G onto Q:
    H = Q @ G @ Q.T
    inside_sq = float(np.square(H).sum())
    outside = max(0.0, total_norm_sq - inside_sq)
    floor = np.sqrt(outside)
    print(f"B={B_keep:3d} -> Subspace projection floor: {floor*1000:7.3f} mHa")
