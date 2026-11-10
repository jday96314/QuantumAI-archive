"""Test Exchange Engine: LOO pruning + Tensor Power replacement on B=719 atoms."""

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

sub = np.load("dfthcsos/submission.npz")
tree = sub["tree_edges"]
words = sub["angle_words"].reshape(-1, 75)
words_2body = words[76:]
U_q = reconstruct_u_vectors(words_2body / 128.0, tree)
all_words = list(words_2body)
all_U = list(torch.tensor(U_q / np.linalg.norm(U_q, axis=1, keepdims=True), dtype=torch.float64, device=device))

mask_diag = ii != jj

def eval_floor(U_list):
    U_curr = torch.stack(U_list, dim=0)
    P = U_curr[:, ii] * U_curr[:, jj] * fac
    S = P @ P.T + 1e-6 * torch.eye(len(U_list), dtype=torch.float64, device=device)
    se, sv = torch.linalg.eigh(S)
    pos = se > 1e-10
    iroot = (sv[:, pos] / torch.sqrt(se[pos])) @ sv[:, pos].T
    Q = iroot @ P
    H = Q @ G @ Q.T
    E = torch.sum(H**2).item()
    floor = np.sqrt(max(0.0, total_norm_sq - E))
    return floor, P, S, Q, H

floor, P, S, Q, H = eval_floor(all_U)
print(f"Initial B={len(all_words)}: Floor = {floor*1000:.4f} mHa")

t0 = time.time()
print("\n--- Running 10 Exchange Cycles ---")
for cycle in range(10):
    B = len(all_U)
    
    # 1. Compute exact LOO loss by evaluating 16 candidate weakest atoms
    # Use Sherman-Morrison / inverse proxy to find candidates:
    Si = torch.linalg.inv(S)
    cap = Si @ H
    K = cap @ Si
    d = Si.diag()
    # Approx loss for each atom
    approx_loss = (2.0 * torch.sum(cap * K.T, dim=1) - (torch.sum(K * cap.T, dim=1))) / d
    candidate_removals = torch.argsort(approx_loss)[:16].cpu().numpy()
    
    # Find the true weakest atom among candidates:
    best_drop_idx = None
    best_drop_E = -1.0
    for cand in candidate_removals:
        keep = [i for i in range(B) if i != cand]
        Uk = [all_U[i] for i in keep]
        fk, _, _, _, _ = eval_floor(Uk)
        if fk < best_drop_E or best_drop_E < 0:
            best_drop_E = fk
            best_drop_idx = cand
            
    # Drop the weakest atom
    dropped_floor = best_drop_E
    keep = [i for i in range(B) if i != best_drop_idx]
    all_words = [all_words[i] for i in keep]
    all_U = [all_U[i] for i in keep]
    floor, P, S, Q, H = eval_floor(all_U)
    
    # 2. Generate a power atom to replace it
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
    best_word = words_cand[best_idx]
    best_u = u_quant_t[best_idx]
    
    all_words.append(best_word)
    all_U.append(best_u)
    floor, P, S, Q, H = eval_floor(all_U)
    
    print(f"Cycle {cycle+1:2d} (B={len(all_words):3d}): swapped atom {best_drop_idx:3d} (drop floor: {dropped_floor*1000:7.4f}) -> New Floor: {floor*1000:7.4f} mHa | Time: {time.time()-t0:.1f}s")
