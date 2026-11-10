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
w = sub["W"][0] # (719, 350)
row_norms = np.linalg.norm(w, axis=1)
sort_order = np.argsort(row_norms) # smallest first

U_q = reconstruct_u_vectors(words_2body / 128.0, tree)
U_t = torch.tensor(U_q / np.linalg.norm(U_q, axis=1, keepdims=True), dtype=torch.float64, device=device)

def get_floor(keep_idx):
    Uk = U_t[keep_idx]
    Pk = Uk[:, ii] * Uk[:, jj] * fac
    Sk = Pk @ Pk.T + 1e-6 * torch.eye(len(keep_idx), dtype=torch.float64, device=device)
    se, sv = torch.linalg.eigh(Sk)
    pos = se > 1e-10
    iroot = (sv[:, pos] / torch.sqrt(se[pos])) @ sv[:, pos].T
    Q = iroot @ Pk
    H = Q @ G @ Q.T
    E = torch.sum(H**2).item()
    return np.sqrt(max(0.0, total_norm_sq - E))

print(f"Full 719 floor: {get_floor(np.arange(719))*1000:.4f} mHa")
for drop_n in [1, 2, 5, 10, 15, 20]:
    keep = sort_order[drop_n:]
    print(f"Drop {drop_n:2d} smallest-norm atoms (B={len(keep):3d}): floor = {get_floor(keep)*1000:7.4f} mHa")
