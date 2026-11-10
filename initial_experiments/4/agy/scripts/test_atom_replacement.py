"""Test replacing lowest-weight atoms in incumbent dictionary with Quartic Tensor Power atoms."""

import time, torch, numpy as np, sys
from pathlib import Path

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
from compress_chem import encode

device = torch.device("cuda")

# 1. Load Reference Spectrum
sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
ii = torch.tensor(sp["ii"], device=device)
jj = torch.tensor(sp["jj"], device=device)
fac = torch.tensor(sp["fac"], dtype=torch.float64, device=device)
ev = torch.tensor(sp["eigenvalues"], dtype=torch.float64, device=device)
V = torch.tensor(sp["eigenvectors"], dtype=torch.float64, device=device)
G = (V * torch.clamp(ev, min=0.0)) @ V.T
total_norm_sq = torch.sum(torch.clamp(ev, min=0.0)**2).item()

# 2. Load incumbent
sub = np.load("dfthcsos/submission.npz")
tree = sub["tree_edges"]
words = sub["angle_words"].reshape(-1, 75)
words_2body = words[76:]
B = len(words_2body) # 719
w_matrix = sub["W"][0] # (719, 350)
weights_norm = np.linalg.norm(w_matrix, axis=1) # (719,)

print(f"Loaded incumbent: B={B}, weights_norm range: [{np.min(weights_norm):.4f}, {np.max(weights_norm):.4f}]")

# Sort atoms by weight norm (lowest first)
sorted_indices = np.argsort(weights_norm)
print(f"Top 5 weakest atom indices: {sorted_indices[:5]} with norms: {weights_norm[sorted_indices[:5]]}")

# Function to compute floor
def compute_floor(word_list):
    U_q = reconstruct_u_vectors(np.array(word_list) / 128.0, tree)
    U_t = torch.tensor(U_q / np.linalg.norm(U_q, axis=1, keepdims=True), dtype=torch.float64, device=device)
    P = U_t[:, ii] * U_t[:, jj] * fac
    S = P @ P.T + 1e-6 * torch.eye(len(word_list), dtype=torch.float64, device=device)
    se, sv = torch.linalg.eigh(S)
    pos = se > 1e-10
    iroot = (sv[:, pos] / torch.sqrt(se[pos])) @ sv[:, pos].T
    Q = iroot @ P
    H = Q @ G @ Q.T
    E = torch.sum(H**2).item()
    fl = np.sqrt(max(0.0, total_norm_sq - E))
    return fl, Q, H

floor_orig, Q, H = compute_floor(words_2body)
print(f"Incumbent floor (B=719): {floor_orig * 1000:.4f} mHa")

# Test dropping k weakest atoms
for k in [1, 5, 10, 20]:
    keep_indices = sorted_indices[k:]
    fl_drop, _, _ = compute_floor([words_2body[i] for i in keep_indices])
    print(f"Dropping {k} weakest atoms (B={719-k}): Floor = {fl_drop * 1000:.4f} mHa (delta: +{(fl_drop - floor_orig)*1000:.4f} mHa)")
