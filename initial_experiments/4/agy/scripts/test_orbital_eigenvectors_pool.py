"""Test if leading orbital-pair eigenvectors of G provide a low projection floor when quantized."""

import numpy as np
import torch
import sys
from pathlib import Path

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
from compress_chem import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

device = torch.device("cuda")
sp = np.load("../../3/parameters/reference_spectrum.npz")
ii = sp["ii"]
jj = sp["jj"]
fac = sp["fac"]
ev = sp["eigenvalues"]
evec = sp["eigenvectors"] # shape (2926, 2926)
N = 76

# Build weighted pool of eigenvectors of V_k
K = 350
pool_vecs = []
pool_weights = []

for k in range(K):
    if ev[k] <= 0:
        continue
    vk = evec[:, k]
    M_k = np.zeros((N, N))
    M_k[ii, jj] = vk / fac
    M_k[jj, ii] = vk / fac
    w, v = np.linalg.eigh(M_k)
    weights = np.sqrt(ev[k]) * np.abs(w)
    for j in range(N):
        if weights[j] > 1e-3:
            pool_vecs.append(v[:, j])
            pool_weights.append(weights[j])

pool_vecs = np.array(pool_vecs)
pool_weights = np.array(pool_weights)
print(f"Pool size: {len(pool_vecs)}")

# Sort pool by weight descending
idx = np.argsort(-pool_weights)
pool_vecs = pool_vecs[idx]
pool_weights = pool_weights[idx]

# Select distinct atoms (filter out near-duplicates)
chosen = []
for v in pool_vecs:
    # check overlap with chosen
    if len(chosen) == 0:
        chosen.append(v)
    else:
        # compute max overlap with chosen
        overlaps = np.abs(np.array(chosen) @ v)
        if np.max(overlaps) < 0.95:
            chosen.append(v)
    if len(chosen) >= 800:
        break

chosen = np.array(chosen)
print(f"Selected {len(chosen)} distinct atoms with mutual overlap < 0.95.")

# Now test quantized projection floor for subsets of chosen:
tree = np.load("../../3/candidates/chem35c45_rot.npz")["tree_edges"]

# Precompute G on GPU
G_gpu = torch.tensor((evec * np.maximum(ev, 0)) @ evec.T, dtype=torch.float64, device=device)
ii_t = torch.tensor(ii, device=device)
jj_t = torch.tensor(jj, device=device)
fac_t = torch.tensor(fac, dtype=torch.float64, device=device)
total_norm_sq = float(np.sum(np.maximum(ev, 0)**2))

for B in [350, 400, 500, 600, 719, len(chosen)]:
    if B > len(chosen):
        continue
    subset = chosen[:B]
    words = encode(subset, tree, 8)
    U_q = reconstruct_u_vectors(words / 128.0, tree)
    U_t = torch.tensor(U_q / np.linalg.norm(U_q, axis=1, keepdims=True), dtype=torch.float64, device=device)
    P = U_t[:, ii_t] * U_t[:, jj_t] * fac_t
    S = P @ P.T
    se, sv = torch.linalg.eigh(S)
    pos = se > 1e-10
    sv = sv[:, pos]
    se = se[pos]
    iroot = (sv / torch.sqrt(se)) @ sv.T
    Q = iroot @ P
    H = Q @ G_gpu @ Q.T
    in_sq = torch.sum(H**2).item()
    floor = np.sqrt(max(0.0, total_norm_sq - in_sq))
    print(f"B = {B:3d} (quantized brot=8) -> floor = {floor*1000:7.3f} mHa")
