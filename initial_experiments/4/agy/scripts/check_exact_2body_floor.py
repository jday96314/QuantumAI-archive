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
U_raw = sub["U"][0]
print(f"sub['U'] shape: {sub['U'].shape}")
print(f"sub['angle_words'] shape: {words.shape}")

# In dfthcsos, words has shape (795, 75): 76 one-body words + 719 two-body words!
words_2body = words[76:]
U_q = reconstruct_u_vectors(words_2body / 128.0, tree)
U_t = torch.tensor(U_q / np.linalg.norm(U_q, axis=1, keepdims=True), dtype=torch.float64, device=device)
B = len(words_2body)

P = U_t[:, ii] * U_t[:, jj] * fac
S = P @ P.T + 1e-6 * torch.eye(B, dtype=torch.float64, device=device)
L = torch.linalg.cholesky(S)
H = P @ G @ P.T
Y = torch.linalg.solve_triangular(L, H, upper=False)
Y = torch.linalg.solve_triangular(L, Y.T, upper=False).T
E_full = torch.sum(Y**2).item()
floor_full = np.sqrt(max(0.0, total_norm_sq - E_full))
print(f"Two-body dictionary B={B}: floor = {floor_full*1000:.4f} mHa")
