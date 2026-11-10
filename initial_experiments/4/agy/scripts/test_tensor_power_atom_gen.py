"""Test Quartic Tensor Power Method on GPU for residual-directed atom generation."""

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

# Let's take an existing dictionary, say the top 400 atoms of submission719
sub = np.load("dfthcsos/submission.npz")
U_sub = torch.tensor(sub["U"][0][:400], dtype=torch.float64, device=device) # (400, 76)
B = 400

P = U_sub[:, ii] * U_sub[:, jj] * fac
S = P @ P.T
se, sv = torch.linalg.eigh(S)
pos = se > 1e-10
iroot = (sv[:, pos] / torch.sqrt(se[pos])) @ sv[:, pos].T
Q = iroot @ P
H = Q @ G @ Q.T
E0 = torch.sum(H**2).item()
floor0 = np.sqrt(max(0.0, total_norm_sq - E0))
print(f"Starting dictionary B={B}: floor = {floor0*1000:.3f} mHa")

# Compute residual projector: Pi_V = Q^T @ Q
# Residual operator: R_perp = G - Q^T @ H @ Q - (I - Q^T Q) G Q^T Q - ...
# More simply: for any vector u, p_perp = (I - Q^T Q) p(u)
# And R_perp @ p = (I - Q^T Q) @ G @ (I - Q^T Q) @ p

# Batch of 128 random starting unit vectors
N_starts = 128
u_batch = torch.randn(N_starts, 76, dtype=torch.float64, device=device)
u_batch = u_batch / torch.norm(u_batch, dim=-1, keepdim=True)

# Precompute mapping from 2926 packed indices to 76x76 matrices
# For each vector v in 2926, unpack(v) is (v / fac)[ii, jj] symmetrized
print("Running power iteration on 128 random seeds...")
t0 = time.time()
for iter in range(10):
    p_batch = u_batch[:, ii] * u_batch[:, jj] * fac # (N_starts, 2926)
    # p_perp = p - (p @ Q.T) @ Q
    q = p_batch @ Q.T
    p_perp = p_batch - q @ Q
    # y = G @ p_perp.T -> (2926, N_starts)
    # y_perp = y - Q.T @ (Q @ y)
    y = G @ p_perp.T
    y_perp = y - Q.T @ (Q @ y) # (2926, N_starts)
    
    # Unpack y_perp for each seed to a 76x76 matrix M, and contract with u:
    # (M u)_p = sum_q M_pq u_q
    # Vectorized:
    # y_sym has shape (N_starts, 76, 76)
    y_scaled = y_perp / fac[:, None] # (2926, N_starts)
    # Compute next u:
    # For each atom m: M_m[p, q] = y_scaled[k, m]
    # (M_m @ u_m)_p = sum_q M_m[p, q] * u_m[q]
    # We can do this with scatter_add!
    # Each pair (ii[k], jj[k]) contributes:
    # p = ii[k]: contributes y_scaled[k, m] * u_m[jj[k]]
    # p = jj[k]: contributes y_scaled[k, m] * u_m[ii[k]] (if ii != jj)
    u_next = torch.zeros_like(u_batch)
    # term 1: to ii from jj
    u_next.index_add_(1, ii, (y_scaled * u_batch[:, jj].T).T)
    # term 2: to jj from ii (where ii != jj)
    mask = ii != jj
    u_next.index_add_(1, jj[mask], (y_scaled[mask] * u_batch[:, ii[mask]].T).T)
    
    # Normalize
    u_batch = u_next / torch.norm(u_next, dim=-1, keepdim=True)

torch.cuda.synchronize()
t1 = time.time()
print(f"10 power iterations on 128 seeds completed in {(t1 - t0)*1000:.1f} ms!")

# Now evaluate exact gains of all 128 generated continuous vectors
A = Q @ G
p_batch = u_batch[:, ii] * u_batch[:, jj] * fac
q = Q @ p_batch.T
norm_p_sq = torch.sum(p_batch**2, dim=1)
norm_q_sq = torch.sum(q**2, dim=0)
d = torch.clamp(norm_p_sq - norm_q_sq, min=1e-10)
Ap = p_batch @ A.T
Hq = H @ q
r = Ap.T - Hq
norm_r_sq = torch.sum(r**2, dim=0)
pGp = torch.sum((p_batch @ G) * p_batch, dim=1)
qAp = torch.sum(q * Ap.T, dim=0)
qHq = torch.sum(q * Hq, dim=0)
s = pGp - 2.0 * qAp + qHq
gains_cont = 2.0 * norm_r_sq / d + (s / d)**2
best_cont_idx = torch.argmax(gains_cont).item()
print(f"Best continuous atom gain: {gains_cont[best_cont_idx].item():.6f}")

# Now quantize ALL 128 generated atoms to brot=8
u_np = u_batch.cpu().numpy()
words = encode(u_np, tree, 8)
u_quant = reconstruct_u_vectors(words / 128.0, tree)
u_quant_t = torch.tensor(u_quant / np.linalg.norm(u_quant, axis=1, keepdims=True), dtype=torch.float64, device=device)

# Evaluate exact gains of quantized atoms!
p_q = u_quant_t[:, ii] * u_quant_t[:, jj] * fac
q_q = Q @ p_q.T
norm_pq_sq = torch.sum(p_q**2, dim=1)
norm_qq_sq = torch.sum(q_q**2, dim=0)
d_q = torch.clamp(norm_pq_sq - norm_qq_sq, min=1e-10)
Ap_q = p_q @ A.T
Hq_q = H @ q_q
r_q = Ap_q.T - Hq_q
norm_rq_sq = torch.sum(r_q**2, dim=0)
pGp_q = torch.sum((p_q @ G) * p_q, dim=1)
qAp_q = torch.sum(q_q * Ap_q.T, dim=0)
qHq_q = torch.sum(q_q * Hq_q, dim=0)
s_q = pGp_q - 2.0 * qAp_q + qHq_q
gains_q = 2.0 * norm_rq_sq / d_q + (s_q / d_q)**2

best_q_idx = torch.argmax(gains_q).item()
print(f"Best quantized atom gain:    {gains_q[best_q_idx].item():.6f}")
print(f"Quantized gain preservation: {gains_q[best_q_idx].item() / gains_cont[best_cont_idx].item() * 100:.1f}%")

# Floor if we add this best quantized atom:
new_floor = np.sqrt(max(0.0, total_norm_sq - (E0 + gains_q[best_q_idx].item())))
print(f"Floor after adding 1 quantized atom: {new_floor*1000:.3f} mHa (drop = {(floor0 - new_floor)*1000:.4f} mHa)")
