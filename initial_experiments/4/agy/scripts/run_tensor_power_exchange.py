"""Tensor Frank-Wolfe with Nonmonotonic Atom Exchange (TFW-AE) on GPU.

Replaces the weakest atoms with newly generated Quartic Tensor Power atoms on the residual,
with periodic discrete STE angle polishing.
"""

import os, sys, time, json
from pathlib import Path

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
import numpy as np
import torch
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

stages = []
for j, (pa, ch) in enumerate(tree):
    pa, ch = int(pa), int(ch)
    if not stages or any(pa in (x, y) or ch in (x, y) for _, x, y in stages[-1]):
        stages.append([])
    stages[-1].append((j, pa, ch))
stages = [tuple(torch.tensor(x, device=device) for x in zip(*s)) for s in stages]

bits = 8
alpha = np.pi / (1 << (bits - 1))
mask_diag = ii != jj

def decode_ste(th):
    angles = (th + (th.round() - th).detach()) * alpha
    un = torch.zeros((len(th), 76), dtype=torch.float64, device=device)
    un[:, 0] = 1.0
    for ts, ps, cs in stages:
        ang = angles[:, ts]
        c, s = ang.cos(), ang.sin()
        up, uc = un[:, ps], un[:, cs]
        new = un.clone()
        new[:, ps] = c * up - s * uc
        new[:, cs] = s * up + c * uc
        un = new
    return un

# Start with ste450_best or top 450 from submission719
init_path = "dfthcsos/ste450_best.npz"
if os.path.exists(init_path):
    z = np.load(init_path)
    words = z["angle_words"]
else:
    sub = np.load("dfthcsos/submission.npz")
    words = sub["angle_words"][76:76+450]

B = len(words)
print(f"Starting TFW-AE with B={B} atoms...")

words_list = [w.copy() for w in words]

def compute_state(words_curr):
    w_arr = np.array(words_curr)
    U_q = reconstruct_u_vectors(w_arr / 128.0, tree)
    U_t = torch.tensor(U_q / np.linalg.norm(U_q, axis=1, keepdims=True), dtype=torch.float64, device=device)
    P = U_t[:, ii] * U_t[:, jj] * fac
    S = P @ P.T + 1e-6 * torch.eye(len(words_curr), dtype=torch.float64, device=device)
    L = torch.linalg.cholesky(S)
    H = P @ G @ P.T
    Y = torch.linalg.solve_triangular(L, H, upper=False)
    Y = torch.linalg.solve_triangular(L, Y.T, upper=False).T
    E = torch.sum(Y**2).item()
    fl = np.sqrt(max(0.0, total_norm_sq - E))
    return U_t, P, S, H, L, E, fl

U_t, P, S, H, L, E, fl = compute_state(words_list)
print(f"Initial discrete floor: {fl*1000:.3f} mHa")

best_floor = fl
best_words = np.array(words_list)
start_time = time.time()

num_swaps = 100
for cycle in range(num_swaps):
    # 1. Identify weakest atom via LOO
    Si = torch.linalg.inv(S)
    cap = Si @ H
    K = cap @ Si
    loss_vectorized = (2.0 * torch.sum(cap * K.T, dim=1) - (K.diag()**2) / Si.diag()) / Si.diag()
    
    # Sort atoms by loss ascending
    weak_indices = torch.argsort(loss_vectorized)
    weak_idx = weak_indices[0].item()
    weak_loss = loss_vectorized[weak_idx].item()
    
    # 2. Temporary remove weak_idx to form basis Q
    keep = [i for i in range(B) if i != weak_idx]
    P_base = P[keep]
    S_base = P_base @ P_base.T + 1e-6 * torch.eye(B - 1, dtype=torch.float64, device=device)
    se, sv = torch.linalg.eigh(S_base)
    pos = se > 1e-10
    iroot = (sv[:, pos] / torch.sqrt(se[pos])) @ sv[:, pos].T
    Q = iroot @ P_base
    H_base = Q @ G @ Q.T
    A = Q @ G
    
    # 3. Quartic Tensor Power Iteration against residual of (B-1) atoms
    N_starts = 64
    u_batch = torch.randn(N_starts, 76, dtype=torch.float64, device=device)
    u_batch = u_batch / torch.norm(u_batch, dim=-1, keepdim=True)
    
    for _ in range(8):
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
        
    # 4. Quantize candidates to 8 bits
    u_np = u_batch.cpu().numpy()
    words_cand = encode(u_np, tree, 8)
    u_quant = reconstruct_u_vectors(words_cand / 128.0, tree)
    u_quant_t = torch.tensor(u_quant / np.linalg.norm(u_quant, axis=1, keepdims=True), dtype=torch.float64, device=device)
    
    # 5. Evaluate exact gains
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
    
    best_cand_idx = torch.argmax(gains_q).item()
    best_gain = gains_q[best_cand_idx].item()
    new_word = words_cand[best_cand_idx]
    
    # Check net improvement
    net_gain = best_gain - weak_loss
    if net_gain > 1e-7:
        # Perform swap!
        words_list[weak_idx] = new_word
        U_t, P, S, H, L, E, fl = compute_state(words_list)
        if fl < best_floor:
            best_floor = fl
            best_words = np.array(words_list)
        status = "SWAP"
    else:
        status = "REJECT"
        
    if (cycle + 1) % 10 == 0:
        elapsed = time.time() - start_time
        print(f"Cycle {cycle+1:3d} | {status} (net = {net_gain:+.6f}) | Floor: {fl*1000:7.3f} mHa (best: {best_floor*1000:7.3f} mHa) | Elapsed: {elapsed:.1f}s")
        
        # Periodic STE Polish every 20 cycles
        if (cycle + 1) % 20 == 0:
            print("  --> Running 15 steps of STE Adam polish...")
            theta_polish = torch.tensor(words_list, dtype=torch.float64, device=device, requires_grad=True)
            opt_p = torch.optim.Adam([theta_polish], lr=0.03)
            for _ in range(15):
                opt_p.zero_grad()
                un_p = decode_ste(theta_polish)
                Pp = un_p[:, ii] * un_p[:, jj] * fac
                Sp = Pp @ Pp.T + 1e-6 * torch.eye(B, dtype=torch.float64, device=device)
                Lp = torch.linalg.cholesky(Sp)
                Hp = Pp @ G @ Pp.T
                Yp = torch.linalg.solve_triangular(Lp, Hp, upper=False)
                Yp = torch.linalg.solve_triangular(Lp, Yp.T, upper=False).T
                (-torch.sum(Yp**2)).backward()
                opt_p.step()
            words_list = [w.copy() for w in theta_polish.round().detach().cpu().numpy().astype(np.int64)]
            U_t, P, S, H, L, E, fl = compute_state(words_list)
            if fl < best_floor:
                best_floor = fl
                best_words = np.array(words_list)
            print(f"  --> Post-polish floor: {fl*1000:7.3f} mHa (best: {best_floor*1000:7.3f} mHa)")

out_path = "/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/dfthcsos/tfw_ae_b450_best.npz"
np.savez(
    out_path,
    angle_words=best_words,
    tree_edges=tree,
    brot=8,
    floor=best_floor,
)
print(f"\nTFW-AE finished in {time.time()-start_time:.1f}s. Saved best candidate to {out_path} with floor = {best_floor*1000:.3f} mHa")
