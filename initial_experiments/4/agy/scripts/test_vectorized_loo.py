"""Test vectorized leave-one-out loss formula on GPU via Sherman-Morrison."""

import time, torch, numpy as np

device = torch.device("cuda")
sp = np.load("../../3/parameters/reference_spectrum.npz")
ii = torch.tensor(sp["ii"], device=device)
jj = torch.tensor(sp["jj"], device=device)
fac = torch.tensor(sp["fac"], dtype=torch.float64, device=device)
ev = torch.tensor(sp["eigenvalues"], dtype=torch.float64, device=device)
V = torch.tensor(sp["eigenvectors"], dtype=torch.float64, device=device)
G = (V * torch.clamp(ev, min=0.0)) @ V.T
total_norm_sq = torch.sum(torch.clamp(ev, min=0.0)**2).item()
tree = np.load("../../3/candidates/chem35c45_rot.npz")["tree_edges"]

# Load ste450_best
z = np.load("dfthcsos/ste450_best.npz")
words = z["angle_words"]
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
from compress_chem import encode
U_q = reconstruct_u_vectors(words / 128.0, tree)
U_t = torch.tensor(U_q / np.linalg.norm(U_q, axis=1, keepdims=True), dtype=torch.float64, device=device)
B = len(words)

P = U_t[:, ii] * U_t[:, jj] * fac
S = P @ P.T + 1e-6 * torch.eye(B, dtype=torch.float64, device=device)
H = P @ G @ P.T
Si = torch.linalg.inv(S)
cap = Si @ H # (B, B)
K = cap @ Si # (B, B)

# Full captured energy
E_full = torch.sum(cap * cap.T).item()
floor_full = np.sqrt(max(0.0, total_norm_sq - E_full))
print(f"Full dictionary B={B}: floor = {floor_full*1000:.3f} mHa")

# Vectorized formula for loss in E when dropping atom b:
# loss_b = (2 * (cap * K.T).sum(1) - (K.diag()**2) / Si.diag()) / Si.diag()
# Let's compute loss for all B atoms:
loss_vectorized = (2.0 * torch.sum(cap * K.T, dim=1) - (K.diag()**2) / Si.diag()) / Si.diag()

# Compare against exact brute-force leave-one-out for 5 random atoms
print("\nVerifying against exact brute-force Cholesky on 5 sample atoms:")
for b in [0, 10, 50, 100, 200]:
    keep = [i for i in range(B) if i != b]
    Pk = P[keep]
    Sk = Pk @ Pk.T + 1e-6 * torch.eye(B - 1, dtype=torch.float64, device=device)
    Lk = torch.linalg.cholesky(Sk)
    Hk = Pk @ G @ Pk.T
    Yk = torch.linalg.solve_triangular(Lk, Hk, upper=False)
    Yk = torch.linalg.solve_triangular(Lk, Yk.T, upper=False).T
    Ek = torch.sum(Yk**2).item()
    actual_loss = E_full - Ek
    pred_loss = loss_vectorized[b].item()
    print(f"  Atom {b:3d}: predicted loss = {pred_loss:.6f}, actual loss = {actual_loss:.6f}, diff = {abs(pred_loss - actual_loss):.2e}")

# Find weakest 20 atoms
weakest = torch.argsort(loss_vectorized)[:20]
print(f"\nWeakest 10 atoms: indices = {weakest[:10].tolist()}")
print(f"Their marginal losses: {loss_vectorized[weakest[:10]].cpu().numpy().round(6)}")
