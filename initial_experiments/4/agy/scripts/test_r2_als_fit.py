"""Solve for optimal block-diagonal (w0, w1) in R=2 on GPU."""

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
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

device = torch.device("cuda")
sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
ii = torch.tensor(sp["ii"], device=device)
jj = torch.tensor(sp["jj"], device=device)
fac = torch.tensor(sp["fac"], dtype=torch.float64, device=device)
ev = torch.tensor(sp["eigenvalues"], dtype=torch.float64, device=device)
V = torch.tensor(sp["eigenvectors"], dtype=torch.float64, device=device)
G = (V * torch.clamp(ev, min=0.0)) @ V.T
total_norm_sq = torch.sum(torch.clamp(ev, min=0.0)**2).item()

sub = np.load("dfthcsos/submission.npz")
tree = sub["tree_edges"]
words = sub["angle_words"].reshape(-1, 75)
words_2body = words[76:]
U_q = reconstruct_u_vectors(words_2body / 128.0, tree)
U_t = torch.tensor(U_q / np.linalg.norm(U_q, axis=1, keepdims=True), dtype=torch.float64, device=device)

# Take 2 ensembles of 350
U0 = U_t[:350]
U1 = U_t[350:700]

P0 = U0[:, ii] * U0[:, jj] * fac # (350, 2926)
P1 = U1[:, ii] * U1[:, jj] * fac # (350, 2926)

# We want to solve for w0 in R^(350 x C) and w1 in R^(350 x C) such that:
# X0 = P0^T w0, X1 = P1^T w1
# M_approx = X0 X0^T + X1 X1^T approx G
# In feature space, P0 G P0^T can be projected.
# Notice that if we split G into two positive parts G0 + G1 = G,
# each G_r can be fitted by X_r X_r^T.
# Let's test the optimal alternating projection / least squares!

# Feature Gram matrices
S0 = P0 @ P0.T + 1e-5 * torch.eye(350, dtype=torch.float64, device=device)
S1 = P1 @ P1.T + 1e-5 * torch.eye(350, dtype=torch.float64, device=device)
L0 = torch.linalg.cholesky(S0)
L1 = torch.linalg.cholesky(S1)

# Alternating least squares:
# Step 0: G0 = 0.5 * G, G1 = 0.5 * G
# Solve w0 to fit G - X1 X1^T, solve w1 to fit G - X0 X0^T
C = 175
X1 = torch.zeros(2926, C, dtype=torch.float64, device=device)

print(f"Solving for optimal R=2 (B=350, C={C}) factor...")
for it in range(30):
    # Target for 0: G_res0 = G - X1 @ X1.T
    G_res0 = G - X1 @ X1.T
    H0 = P0 @ G_res0 @ P0.T
    h0_sym = torch.linalg.solve_triangular(L0, H0, upper=False)
    h0_sym = torch.linalg.solve_triangular(L0, h0_sym.T, upper=False).T
    e0, v0 = torch.linalg.eigh(0.5 * (h0_sym + h0_sym.T))
    top_e0 = torch.clamp(e0[-C:], min=0.0)
    top_v0 = v0[:, -C:]
    w0 = torch.linalg.solve_triangular(L0.T, top_v0 * torch.sqrt(top_e0), upper=True)
    X0 = P0.T @ w0
    
    # Target for 1: G_res1 = G - X0 @ X0.T
    G_res1 = G - X0 @ X0.T
    H1 = P1 @ G_res1 @ P1.T
    h1_sym = torch.linalg.solve_triangular(L1, H1, upper=False)
    h1_sym = torch.linalg.solve_triangular(L1, h1_sym.T, upper=False).T
    e1, v1 = torch.linalg.eigh(0.5 * (h1_sym + h1_sym.T))
    top_e1 = torch.clamp(e1[-C:], min=0.0)
    top_v1 = v1[:, -C:]
    w1 = torch.linalg.solve_triangular(L1.T, top_v1 * torch.sqrt(top_e1), upper=True)
    X1 = P1.T @ w1
    
    # Compute error:
    M_approx = X0 @ X0.T + X1 @ X1.T
    diff = M_approx - G
    err2 = torch.sqrt(torch.sum(diff * diff)).item()
    if (it + 1) % 5 == 0 or it == 0:
        print(f"Iter {it+1:2d}: err2 = {err2*1000:7.3f} mHa")
