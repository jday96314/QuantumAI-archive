"""Test growing an 8-bit quantized dictionary from scratch using Quartic Tensor Power."""

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
from compress_chem import encode

device = torch.device("cuda")
sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
ii = torch.tensor(sp["ii"], device=device)
jj = torch.tensor(sp["jj"], device=device)
fac = torch.tensor(sp["fac"], dtype=torch.float64, device=device)
ev = torch.tensor(sp["eigenvalues"], dtype=torch.float64, device=device)
V = torch.tensor(sp["eigenvectors"], dtype=torch.float64, device=device)
G = (V * torch.clamp(ev, min=0.0)) @ V.T
total_norm_sq = torch.sum(torch.clamp(ev, min=0.0)**2).item()
tree = np.load("../../3/candidates/chem35c45_rot.npz")["tree_edges"]

mask_diag = ii != jj

# Start with top eigenvectors or seed atoms
# Let's seed with top 50 1-body density vectors or random Givens atoms
print(f"Total reference energy norm: {np.sqrt(total_norm_sq)*1000:.3f} mHa")

# Initial seed: 32 random 8-bit atoms
u_init = torch.randn(32, 76, dtype=torch.float64, device=device)
u_init = (u_init / torch.norm(u_init, dim=-1, keepdim=True)).cpu().numpy()
words_init = encode(u_init, tree, 8)
u_quant = reconstruct_u_vectors(words_init / 128.0, tree)
u_quant_t = torch.tensor(u_quant / np.linalg.norm(u_quant, axis=1, keepdims=True), dtype=torch.float64, device=device)

all_words = list(words_init)
all_U = list(u_quant_t)

U_curr = torch.stack(all_U, dim=0)
P = U_curr[:, ii] * U_curr[:, jj] * fac
S = P @ P.T + 1e-6 * torch.eye(len(all_U), dtype=torch.float64, device=device)
se, sv = torch.linalg.eigh(S)
pos = se > 1e-10
iroot = (sv[:, pos] / torch.sqrt(se[pos])) @ sv[:, pos].T
Q = iroot @ P
H = Q @ G @ Q.T
E = torch.sum(H**2).item()
floor = np.sqrt(max(0.0, total_norm_sq - E))
print(f"Starting seed B={len(all_words)}: Floor = {floor*1000:.3f} mHa")

t0 = time.time()
for step in range(100):
    A = Q @ G
    N_starts = 64
    u_batch = torch.randn(N_starts, 76, dtype=torch.float64, device=device)
    u_batch = u_batch / torch.norm(u_batch, dim=-1, keepdim=True)
    
    for it in range(8):
        p_b = u_batch[:, ii] * u_batch[:, jj] * fac
        q = p_b @ Q.T
        p_perp = p_b - q @ Q
        y = G @ p_perp.T
        y_perp = y - Q.T @ (Q @ y)
        y_scaled = y_perp / fac[:, None]
        
        u_next = torch.zeros_like(u_batch)
        u_next.index_add_(1, ii, (y_scaled * u_batch[:, jj].T).T)
        u_next.index_add_(1, jj[mask_diag], (y_scaled[mask_diag] * u_batch[:, ii[mask_diag]].T).T)
        u_batch = u_next / torch.norm(u_next, dim=-1, keepdim=True)
        
    u_np = u_batch.cpu().numpy()
    words_cand = encode(u_np, tree, 8)
    u_quant = reconstruct_u_vectors(words_cand / 128.0, tree)
    u_quant_t = torch.tensor(u_quant / np.linalg.norm(u_quant, axis=1, keepdims=True), dtype=torch.float64, device=device)
    
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
    
    best_idx = torch.argmax(gains_q).item()
    best_gain = gains_q[best_idx].item()
    best_word = words_cand[best_idx]
    best_u = u_quant_t[best_idx]
    
    all_words.append(best_word)
    all_U.append(best_u)
    
    U_curr = torch.stack(all_U, dim=0)
    P = U_curr[:, ii] * U_curr[:, jj] * fac
    S = P @ P.T + 1e-6 * torch.eye(len(all_U), dtype=torch.float64, device=device)
    se, sv = torch.linalg.eigh(S)
    pos = se > 1e-10
    iroot = (sv[:, pos] / torch.sqrt(se[pos])) @ sv[:, pos].T
    Q = iroot @ P
    H = Q @ G @ Q.T
    E = torch.sum(H**2).item()
    floor = np.sqrt(max(0.0, total_norm_sq - E))
    
    if (step + 1) % 10 == 0:
        print(f"B={len(all_words):3d} | Floor: {floor*1000:7.3f} mHa (gain: {best_gain:.5e}) | Elapsed: {time.time()-t0:.1f}s")
