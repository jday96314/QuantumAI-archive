"""Test rank C projection error for submission719 across C in [300..350]."""

import time, torch, numpy as np
import sys
from pathlib import Path

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]

device = torch.device("cuda")
sp = np.load("../../3/parameters/reference_spectrum.npz")
ii = torch.tensor(sp["ii"], device=device)
jj = torch.tensor(sp["jj"], device=device)
fac = torch.tensor(sp["fac"], dtype=torch.float64, device=device)
ev = torch.tensor(sp["eigenvalues"], dtype=torch.float64, device=device)
V = torch.tensor(sp["eigenvectors"], dtype=torch.float64, device=device)
G = (V * torch.clamp(ev, min=0.0)) @ V.T
total_norm_sq = torch.sum(torch.clamp(ev, min=0.0)**2).item()

sub = np.load("dfthcsos/submission.npz")
U_t = torch.tensor(sub["U"][0], dtype=torch.float64, device=device)
B = len(U_t)

P = U_t[:, ii] * U_t[:, jj] * fac
S = P @ P.T + 1e-6 * torch.eye(B, dtype=torch.float64, device=device)
H = P @ G @ P.T
L = torch.linalg.cholesky(S)
hh = torch.linalg.solve_triangular(L, H, upper=False)
hh = torch.linalg.solve_triangular(L, hh.T, upper=False).T
e, v = torch.linalg.eigh(0.5 * (hh + hh.T))

print(f"Optimal low-rank truncation of submission719:")
for C in [350, 345, 340, 335, 330, 325, 320, 300]:
    # Truncate to top C components
    e_top = e[-C:].clamp(min=0.0)
    v_top = v[:, -C:]
    w = torch.linalg.solve_triangular(L.T, v_top * torch.sqrt(e_top), upper=True) # (B, C)
    K = w @ w.T
    KS = K @ S
    err2 = float(torch.sqrt(torch.clamp(total_norm_sq - 2.0 * torch.sum(K * H) + torch.sum(KS * KS.T), min=0.0)))
    print(f"  C = {C:3d} -> err2 = {err2*1000:7.3f} mHa")
