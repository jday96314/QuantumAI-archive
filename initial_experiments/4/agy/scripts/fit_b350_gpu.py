"""Initialize U from dominant eigenvectors of top H2 slices, then optimize on GPU."""

import sys
import time
from pathlib import Path
import numpy as np
import torch

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
ev, V = sp["eigenvalues"], sp["eigenvectors"]
ii, jj = sp["ii"], sp["jj"]
fac = sp["fac"]

B = 350
N = 76

# Construct top matrices M_k
order = np.argsort(ev)[::-1]
vectors = []
for k in range(B // 4 + 10):
    idx = order[k]
    val = ev[idx]
    if val <= 0:
        break
    vec = V[:, idx] * np.sqrt(val)
    M = np.zeros((N, N))
    M[ii, jj] = vec * fac
    M[jj, ii] = vec * fac
    e, v = np.linalg.eigh(M)
    # Take top eigenvectors of M
    for rank_idx in np.argsort(np.abs(e))[::-1][:4]:
        u_cand = v[:, rank_idx]
        vectors.append(u_cand)
        if len(vectors) >= B:
            break
    if len(vectors) >= B:
        break

U0 = np.array(vectors[:B])
U0 /= np.linalg.norm(U0, axis=1, keepdims=True)
print(f"Initialized U0 with shape {U0.shape}")

# Now optimize U on GPU using PyTorch
device = "cuda"
dtype = torch.float64

u = torch.tensor(U0, device=device, dtype=dtype, requires_grad=True)
G = torch.tensor((V * np.maximum(ev, 0)) @ V.T, device=device, dtype=dtype)
ti = torch.tensor(ii, device=device)
tj = torch.tensor(jj, device=device)
ff = torch.tensor(fac, device=device, dtype=dtype)
eye = torch.eye(B, device=device, dtype=dtype)
const = G.square().sum()

opt = torch.optim.Adam([u], lr=0.003)
best_err = 1e9

start = time.time()
print("Starting Adam optimization of B=350 on RTX 4090...")

for step in range(3000):
    opt.zero_grad()
    un = u / u.norm(dim=1, keepdim=True)
    P = un[:, ti] * un[:, tj] * ff
    S = (un @ un.T).square()
    H = (P @ G) @ P.T
    Sr = S + 1e-6 * eye
    L = torch.linalg.cholesky(Sr)
    K = torch.cholesky_solve(torch.cholesky_solve(H, L).T, L).T
    K = 0.5 * (K + K.T)
    KS = K @ Sr
    loss = const - 2 * (K * H).sum() + (KS * KS.T).sum()
    loss.backward()
    opt.step()

    # Decay lr
    decay = 0.1 + 0.9 * 0.5 * (1 + np.cos(np.pi * step / 3000))
    opt.param_groups[0]["lr"] = 0.003 * decay

    if step % 200 == 0:
        with torch.no_grad():
            physical = K @ S
            err2 = float((const - 2 * (K * H).sum() + (physical * physical.T).sum()).clamp(min=0).sqrt())
            lam_bound = float(0.25 * K.abs().sum())
            if err2 < best_err:
                best_err = err2
            print(f"Step {step:4d} | err2: {err2*1000:7.3f} mHa | best: {best_err*1000:7.3f} mHa | lambda_bound: {lam_bound:.1f} | time: {time.time()-start:.1f}s")

out_file = "/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/runs/fit_b350.npz"
np.savez(out_file, U=un.detach().cpu().numpy()[None, :, :], K=K.detach().cpu().numpy(), err2=best_err)
print(f"Done! Final best err2: {best_err*1000:.3f} mHa. Saved to {out_file}")
