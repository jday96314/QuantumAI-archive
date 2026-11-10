"""Replace k weakest atoms with k tensor power atoms, then evaluate floor and solve weights."""

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
mask_diag = ii != jj

# 2. Load incumbent
sub = np.load("dfthcsos/submission.npz")
tree = sub["tree_edges"]
words = sub["angle_words"].reshape(-1, 75)
words_1body = words[:76]
words_2body = words[76:]
B = len(words_2body) # 719
w_matrix = sub["W"][0] # (719, 350)
weights_norm = np.linalg.norm(w_matrix, axis=1) # (719,)

# Sort atoms by weight norm (lowest first)
sorted_indices = np.argsort(weights_norm)

def power_iterations(current_words, n_new=5, n_starts=128, n_inner=10):
    U_q = reconstruct_u_vectors(np.array(current_words) / 128.0, tree)
    all_U = [torch.tensor(u / np.linalg.norm(u), dtype=torch.float64, device=device) for u in U_q]
    new_words = list(current_words)
    
    t0 = time.time()
    for step in range(n_new):
        U_curr = torch.stack(all_U, dim=0)
        P = U_curr[:, ii] * U_curr[:, jj] * fac
        S = P @ P.T + 1e-6 * torch.eye(len(all_U), dtype=torch.float64, device=device)
        se, sv = torch.linalg.eigh(S)
        pos = se > 1e-10
        iroot = (sv[:, pos] / torch.sqrt(se[pos])) @ sv[:, pos].T
        Q = iroot @ P
        H = Q @ G @ Q.T
        E = torch.sum(H**2).item()
        cur_floor = np.sqrt(max(0.0, total_norm_sq - E))
        
        A = Q @ G
        u_batch = torch.randn(n_starts, 76, dtype=torch.float64, device=device)
        u_batch = u_batch / torch.norm(u_batch, dim=-1, keepdim=True)
        
        for it in range(n_inner):
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
        
        new_words.append(best_word)
        all_U.append(best_u)
        print(f"  Step {step+1}/{n_new}: floor = {cur_floor*1000:.4f} mHa -> gain = {gains_q[best_idx].item()*1e6:.2f} uHa")
        
    U_curr = torch.stack(all_U, dim=0)
    P = U_curr[:, ii] * U_curr[:, jj] * fac
    S = P @ P.T + 1e-6 * torch.eye(len(all_U), dtype=torch.float64, device=device)
    se, sv = torch.linalg.eigh(S)
    pos = se > 1e-10
    iroot = (sv[:, pos] / torch.sqrt(se[pos])) @ sv[:, pos].T
    Q = iroot @ P
    H = Q @ G @ Q.T
    E = torch.sum(H**2).item()
    final_floor = np.sqrt(max(0.0, total_norm_sq - E))
    print(f"Done in {time.time()-t0:.2f}s: Final Floor = {final_floor*1000:.4f} mHa")
    return new_words, final_floor

print("--- Testing Replacing 5 Weakest Atoms with 5 Power Atoms ---")
keep_indices = sorted_indices[5:]
trimmed_words = [words_2body[i] for i in keep_indices]
replaced_words, fl = power_iterations(trimmed_words, n_new=5, n_starts=128, n_inner=10)
print(f"Original B=719 floor: 19.4272 mHa | New B=719 floor: {fl*1000:.4f} mHa")
