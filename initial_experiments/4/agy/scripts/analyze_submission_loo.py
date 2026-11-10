import time, torch, numpy as np
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
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
U_q = reconstruct_u_vectors(words / 128.0, tree)
U_t = torch.tensor(U_q / np.linalg.norm(U_q, axis=1, keepdims=True), dtype=torch.float64, device=device)
B = len(words)

P = U_t[:, ii] * U_t[:, jj] * fac
S = P @ P.T + 1e-6 * torch.eye(B, dtype=torch.float64, device=device)
L = torch.linalg.cholesky(S)
H = P @ G @ P.T
Y = torch.linalg.solve_triangular(L, H, upper=False)
Y = torch.linalg.solve_triangular(L, Y.T, upper=False).T
E_full = torch.sum(Y**2).item()
floor_full = np.sqrt(max(0.0, total_norm_sq - E_full))
print(f"Full dictionary B={B}: floor = {floor_full*1000:.4f} mHa")

# Sherman-Morrison LOO formula
Si = torch.linalg.inv(S)
cap = Si @ H
K = cap @ Si

# Let's verify accuracy of LOO loss
# Exact LOO delta for dropping atom b
drop_losses = []
for b in range(B):
    # Sherman-Morrison exact update:
    # (S_{-b})^{-1} = S_sub^{-1}
    # Or directly from Cholesky solve:
    pass

# Vectorized LOO energy delta
d = Si.diag()
# E_drop delta approx:
delta_E = (2.0 * torch.sum(cap * K.T, dim=1) - (torch.sum(K * cap.T, dim=1))) / d
# Sort atoms by delta_E
sorted_indices = torch.argsort(delta_E)

print("Top 10 least important atoms (smallest delta E):")
for i in range(10):
    idx = sorted_indices[i].item()
    print(f"  Atom {idx:3d}: delta E = {delta_E[idx].item():.8e}")

print("\nTop 10 most important atoms:")
for i in range(1, 11):
    idx = sorted_indices[-i].item()
    print(f"  Atom {idx:3d}: delta E = {delta_E[idx].item():.8e}")

# What is the floor if we drop the 19 least important atoms? (719 -> 700)
for drop_count in [5, 10, 19, 25, 50, 70]:
    keep_indices = sorted_indices[drop_count:].cpu().numpy()
    Pk = P[keep_indices]
    Sk = Pk @ Pk.T + 1e-6 * torch.eye(len(keep_indices), dtype=torch.float64, device=device)
    Lk = torch.linalg.cholesky(Sk)
    Hk = Pk @ G @ Pk.T
    Yk = torch.linalg.solve_triangular(Lk, Hk, upper=False)
    Yk = torch.linalg.solve_triangular(Lk, Yk.T, upper=False).T
    Ek = torch.sum(Yk**2).item()
    floor_k = np.sqrt(max(0.0, total_norm_sq - Ek))
    print(f"Drop {drop_count:2d} atoms (B={len(keep_indices):3d}): floor = {floor_k*1000:.4f} mHa")
