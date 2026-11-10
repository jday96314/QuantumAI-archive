"""Check projection floor of R=2 partition of incumbent 719 atoms: U0=atoms[:350], U1=atoms[350:700]."""

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

# Partition into 2 ensembles of 350
U0 = U_t[:350]
U1 = U_t[350:700]

P0 = U0[:, ii] * U0[:, jj] * fac # (350, 2926)
P1 = U1[:, ii] * U1[:, jj] * fac # (350, 2926)
P_all = torch.cat([P0, P1], dim=0) # (700, 2926)

S = P_all @ P_all.T + 1e-6 * torch.eye(700, dtype=torch.float64, device=device)
se, sv = torch.linalg.eigh(S)
pos = se > 1e-10
iroot = (sv[:, pos] / torch.sqrt(se[pos])) @ sv[:, pos].T
Q = iroot @ P_all
H = Q @ G @ Q.T
E = torch.sum(H**2).item()
floor = np.sqrt(max(0.0, total_norm_sq - E))
print(f"R=2 (B=350, 350) projection floor: {floor*1000:.4f} mHa")
