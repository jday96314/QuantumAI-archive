"""Solve for optimal w given quantized U on chem4b256c128."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
import numpy as np
from compress_chem import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

z = np.load(ROOT / "initial_experiments/3/candidates/chem4b256c128_float_rot_off.npz")
U = z["U"]
w_orig = z["w"]
R, B, N = U.shape
C = w_orig.shape[2]

tree = np.load(ROOT / "initial_experiments/3/candidates/chem35c45_rot.npz")["tree_edges"]

sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
ii, jj, fac = sp["ii"], sp["jj"], sp["fac"]
ev, V = sp["eigenvalues"], sp["eigenvectors"]
G = (V * np.maximum(ev, 0)) @ V.T

print(f"R={R}, B={B}, C={C}, N={N}")

for brot in [8, 9, 10, 11, 12]:
    # Encode U to brot bits
    words = encode(U.reshape(-1, N), tree, brot)
    U_q = reconstruct_u_vectors(words / (1 << (brot - 1)), tree)
    U_reshaped = U_q.reshape(R, B, N)
    u_norm = U_reshaped / np.linalg.norm(U_reshaped, axis=2, keepdims=True)
    
    # Feature matrix for each ensemble r: P_r has shape (B, 2926)
    # Total feature matrix P has shape (R * B, 2926)
    # The represented tensor is sum_r P_r.T @ (w_r @ w_r.T) @ P_r
    # Or in factored form: for each r, w_r has shape (B, C)
    # X_r = P_r.T @ w_r has shape (2926, C)
    # sum_r X_r @ X_r.T = G
    # Let's project G onto the span of P
    P_all = (u_norm[:, :, ii] * u_norm[:, :, jj] * fac).reshape(R * B, -1)
    
    # Projection floor of this quantized dictionary
    S = P_all @ P_all.T
    se, sv = np.linalg.eigh(S)
    pos = se > 1e-10
    sv = sv[:, pos]
    se = se[pos]
    iroot = (sv / np.sqrt(se)) @ sv.T
    Q = iroot @ P_all
    H = Q @ G @ Q.T
    inside_sq = float(np.square(H).sum())
    total_norm_sq = float(np.square(np.maximum(ev, 0)).sum())
    outside = max(0.0, total_norm_sq - inside_sq)
    floor = np.sqrt(outside)
    print(f"brot={brot:2d} -> Quantized dictionary projection floor: {floor*1000:7.3f} mHa")
