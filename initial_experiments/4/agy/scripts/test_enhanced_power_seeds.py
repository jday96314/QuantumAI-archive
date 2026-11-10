"""Test enhanced seed pool for Quartic Tensor Power method on GPU."""

import time, torch, numpy as np
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
ii = torch.tensor(sp["ii"], device=device)
jj = torch.tensor(sp["jj"], device=device)
fac = torch.tensor(sp["fac"], dtype=torch.float64, device=device)
ev = torch.tensor(sp["eigenvalues"], dtype=torch.float64, device=device)
V = torch.tensor(sp["eigenvectors"], dtype=torch.float64, device=device)
G = (V * torch.clamp(ev, min=0.0)) @ V.T
total_norm_sq = torch.sum(torch.clamp(ev, min=0.0)**2).item()
tree = np.load("../../3/candidates/chem35c45_rot.npz")["tree_edges"]

# Load grow650
z = np.load("../../3/candidates/grow650_q8_c350.npz")
words = z["angle_words"].reshape(-1, 75)
B = len(words)
U_q = reconstruct_u_vectors(words / 128.0, tree)
U_t = torch.tensor(U_q / np.linalg.norm(U_q, axis=1, keepdims=True), dtype=torch.float64, device=device)

# Find weakest atom
P = U_t[:, ii] * U_t[:, jj] * fac
S = P @ P.T + 1e-6 * torch.eye(B, dtype=torch.float64, device=device)
H = P @ G @ P.T
Si = torch.linalg.inv(S)
cap = Si @ H
K = cap @ Si
loss_vectorized = (2.0 * torch.sum(cap * K.T, dim=1) - (K.diag()**2) / Si.diag()) / Si.diag()

weak_idx = torch.argmin(loss_vectorized).item()
weak_loss = loss_vectorized[weak_idx].item()
print(f"Weakest atom {weak_idx}: marginal loss = {weak_loss:.6f}")

# Subspace without weak atom
keep = [i for i in range(B) if i != weak_idx]
P_base = P[keep]
S_base = P_base @ P_base.T + 1e-6 * torch.eye(B - 1, dtype=torch.float64, device=device)
se, sv = torch.linalg.eigh(S_base)
pos = se > 1e-10
iroot = (sv[:, pos] / torch.sqrt(se[pos])) @ sv[:, pos].T
Q = iroot @ P_base
H_base = Q @ G @ Q.T
A = Q @ G

# Compute leading eigenvectors of residual tensor R_perp
# R_perp = (I - Q^T Q) G (I - Q^T Q)
QG = Q @ G
R_perp = G - Q.T @ QG - QG.T @ Q + Q.T @ (QG @ Q.T) @ Q
# Top 5 eigenvectors of R_perp
ev_r, v_r = torch.linalg.eigh(R_perp)
print(f"Top 3 eigenvalues of R_perp: {ev_r[-3:].cpu().numpy().round(4)}")

# Extract orbital eigenvectors from top 3 modes of R_perp
seed_vecs = []
# 1. Weak atom itself
seed_vecs.append(U_t[weak_idx])

# 2. Top 3 modes of R_perp
for mode_idx in [-1, -2, -3]:
    top_v = v_r[:, mode_idx]
    M_mode = torch.zeros((76, 76), dtype=torch.float64, device=device)
    M_mode[ii, jj] = top_v / fac
    M_mode[jj, ii] = top_v / fac
    w_m, v_m = torch.linalg.eigh(M_mode)
    # add top 5 positive and top 5 negative eigenvectors
    for j in list(range(5)) + list(range(71, 76)):
        seed_vecs.append(v_m[:, j])

# 3. Random seeds
N_random = 64
u_rand = torch.randn(N_random, 76, dtype=torch.float64, device=device)
for j in range(N_random):
    seed_vecs.append(u_rand[j] / torch.norm(u_rand[j]))

u_seeds = torch.stack(seed_vecs, dim=0)
u_seeds = u_seeds / torch.norm(u_seeds, dim=-1, keepdim=True)
print(f"Total seeds: {len(u_seeds)}")

# Run 20 power iterations
mask_diag = ii != jj
u_curr = u_seeds.clone()
t0 = time.time()
for it in range(20):
    p_b = u_curr[:, ii] * u_curr[:, jj] * fac
    q = p_b @ Q.T
    p_perp = p_b - q @ Q
    y = G @ p_perp.T
    y_perp = y - Q.T @ (Q @ y)
    y_scaled = y_perp / fac[:, None]
    
    u_next = torch.zeros_like(u_curr)
    u_next.index_add_(1, ii, (y_scaled * u_curr[:, jj].T).T)
    u_next.index_add_(1, jj[mask_diag], (y_scaled[mask_diag] * u_curr[:, ii[mask_diag]].T).T)
    u_curr = u_next / torch.norm(u_next, dim=-1, keepdim=True)

torch.cuda.synchronize()
print(f"20 power iterations completed in {(time.time()-t0)*1000:.1f} ms!")

# Quantize and score
u_np = u_curr.cpu().numpy()
words_cand = encode(u_np, tree, 8)
u_quant = reconstruct_u_vectors(words_cand / 128.0, tree)
u_quant_t = torch.tensor(u_quant / np.linalg.norm(u_quant, axis=1, keepdims=True), dtype=torch.float64, device=device)

p_q = u_quant_t[:, ii] * u_quant_t[:, jj] * fac
q_q = Q @ p_q.T
norm_pq_sq = torch.sum(p_q**2, dim=1)
norm_qq_sq = torch.sum(q_q**2, dim=0)
d_q = torch.clamp(norm_pq_sq - norm_qq_sq, min=1e-10)
Ap_q = p_q @ A.T
Hq_q = H_base @ q_q
r_q = Ap_q.T - Hq_q
norm_rq_sq = torch.sum(r_q**2, dim=0)
pGp_q = torch.sum((p_q @ G) * p_q, dim=1)
qAp_q = torch.sum(q_q * Ap_q.T, dim=0)
qHq_q = torch.sum(q_q * Hq_q, dim=0)
s_q = pGp_q - 2.0 * qAp_q + qHq_q
gains_q = 2.0 * norm_rq_sq / d_q + (s_q / d_q)**2

best_gain = torch.max(gains_q).item()
best_idx = torch.argmax(gains_q).item()
print(f"Best quantized gain: {best_gain:.6f}")
print(f"Weak atom loss:     {weak_loss:.6f}")
print(f"Net gain:           {best_gain - weak_loss:+.6f}")
if best_gain > weak_loss:
    print(f"*** SUCCESSFUL SWAP FOUND! Net improvement = {best_gain - weak_loss:+.6f} ***")
else:
    print("Swap rejected.")
