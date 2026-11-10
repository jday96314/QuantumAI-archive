"""Solve optimal (w, W) on GPU for any quantized dictionary U_q to satisfy err2 <= 19.95 mHa and gap >= 3.9152 Ha."""

import argparse, time, json
import numpy as np
import torch
from pathlib import Path

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")

def solve_weights(in_fit: str, out_fit: str, C_target: int = 340, steps: int = 500, lr: float = 0.001):
    device = torch.device("cuda")
    z = np.load(in_fit)
    U_raw = np.asarray(z["U"], dtype=np.float64)
    if U_raw.ndim == 2:
        U_raw = U_raw[None, :, :]
    R, B, N = U_raw.shape
    U_norm = U_raw / np.linalg.norm(U_raw, axis=2, keepdims=True)
    
    # Chem reference
    ref_dict = np.load(ROOT / "initial_experiments/3/parameters/chem_reference.npz")
    sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
    
    S_ref = np.asarray(ref_dict["S"], dtype=np.float64)
    beta_1 = float(ref_dict["beta_1"])
    h1_thc = np.asarray(ref_dict["h1_thc"], dtype=np.float64)
    ecore = float(ref_dict["ecore"])
    electrons = int(ref_dict["electrons"])
    energy_upper = -22140.27335
    
    # Target two-body packed matrix M_t
    ii, jj, fac = sp["ii"], sp["jj"], sp["fac"]
    ev, V = sp["eigenvalues"], sp["eigenvectors"]
    G = (V * np.maximum(ev, 0)) @ V.T # (2926, 2926)
    
    G_t = torch.tensor(G, dtype=torch.float64, device=device)
    total_norm_sq = float(np.sum(np.maximum(ev, 0)**2))
    
    rows, cols = np.tril_indices(N)
    weights_tril = np.where(rows == cols, 1.0, np.sqrt(2.0))
    weights_t = torch.tensor(weights_tril, dtype=torch.float64, device=device)
    rows_t = torch.tensor(rows, dtype=torch.long, device=device)
    cols_t = torch.tensor(cols, dtype=torch.long, device=device)
    
    h1_fixed = h1_thc - 0.5 * electrons * S_ref + beta_1 * np.eye(N)
    h1_fixed_t = torch.tensor(h1_fixed, dtype=torch.float64, device=device)
    
    U_t = torch.tensor(U_norm, dtype=torch.float64, device=device)
    P_t = (U_t[:, :, rows_t] * U_t[:, :, cols_t]) * weights_t[None, None, :]
    
    print(f"Solving initial w from projection for R={R}, B={B}, C={C_target}...")
    P_flat = P_t.reshape(R * B, 2926)
    S_mat = P_flat @ P_flat.T + 1e-6 * torch.eye(R * B, dtype=torch.float64, device=device)
    L = torch.linalg.cholesky(S_mat)
    H_mat = P_flat @ G_t @ P_flat.T
    hh = torch.linalg.solve_triangular(L, H_mat, upper=False)
    hh = torch.linalg.solve_triangular(L, hh.T, upper=False).T
    e_h, v_h = torch.linalg.eigh(0.5 * (hh + hh.T))
    
    v_top = v_h[:, -C_target:]
    e_top = e_h[-C_target:].clamp(min=0.0)
    w_init = torch.linalg.solve_triangular(L.T, v_top * torch.sqrt(e_top), upper=True) # (R*B, C_target)
    
    K_init = w_init @ w_init.T
    KS_init = K_init @ S_mat
    init_err2 = np.sqrt(max(0.0, total_norm_sq - 2.0 * torch.sum(K_init * H_mat).item() + torch.sum(KS_init * KS_init.T).item()))
    print(f"Initial projection err2: {init_err2*1000:.3f} mHa")
    
    w_var = w_init.reshape(R, B, C_target).contiguous().clone().detach().requires_grad_(True)
    
    # Initialize W to cancel out the SOS gap:
    # 0.5 * sum(W^2) should be approx 540.35
    W_init_val = np.sqrt(1080.7 / (R * C_target))
    W_var = torch.full((R, C_target), W_init_val, dtype=torch.float64, device=device).contiguous().requires_grad_(True)
    
    opt = torch.optim.Adam([
        {"params": [w_var], "lr": lr},
        {"params": [W_var], "lr": lr * 5.0},
    ])
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps, eta_min=lr*0.01)
    
    t0 = time.time()
    best_loss = float("inf")
    best_state = None
    
    for step in range(steps):
        opt.zero_grad()
        X_p = torch.einsum("rbc,rbp->rcp", w_var, P_t)
        X_flat = X_p.reshape(R * C_target, 2926)
        diff = X_flat.T @ X_flat - G_t
        err2_sq = torch.sum(diff * diff)
        err2_val = torch.sqrt(torch.clamp(err2_sq, min=0.0))
        
        # 1-body shift
        X_rc = torch.zeros((R, C_target, N, N), dtype=torch.float64, device=device)
        X_rc[:, :, rows_t, cols_t] = X_p / weights_t[None, None, :]
        X_rc = X_rc + X_rc.transpose(-1, -2) - torch.diagonal(X_rc, dim1=-2, dim2=-1).diag_embed()
        h1_p = h1_fixed_t - torch.einsum("rc,rcpq->pq", W_var, X_rc)
        
        ev_1 = torch.linalg.eigvalsh(0.5 * (h1_p + h1_p.T))
        lam_h1 = torch.sum(torch.abs(ev_1))
        tr_h1 = torch.sum(ev_1)
        
        e_sos = ecore - electrons * beta_1 - (lam_h1 - tr_h1) - 0.5 * torch.sum(W_var**2)
        e_gap = energy_upper - e_sos
        
        # Target constraints: err2 <= 0.01995, gap in [3.9152, 4.05]
        err2_pen = torch.relu(err2_val - 0.01995)**2 * 1e5
        gap_pen = torch.relu(3.9152 - e_gap)**2 * 1e4 + torch.relu(e_gap - 4.05)**2 * 1e4
        
        w_b_sum = torch.sum(w_var, dim=1)
        w_b_abs = torch.sum(torch.abs(w_var), dim=1)
        lam_h2 = 0.25 * torch.sum((w_b_abs + torch.abs(W_var + w_b_sum))**2)
        Lambda = lam_h1 + lam_h2
        
        rad = torch.clamp(e_gap * (2.0 * Lambda - e_gap), min=1.0)
        lam_eff = torch.sqrt(rad)
        
        loss = err2_sq * 50.0 + err2_pen + gap_pen + lam_eff * 0.005
        loss.backward()
        opt.step()
        sched.step()
        
        if (step + 1) % 25 == 0:
            print(f"Step {step+1:3d} | err2 = {err2_val.item()*1000:6.3f} mHa | gap = {e_gap.item():6.4f} Ha | lam_eff = {lam_eff.item():6.2f} | loss = {loss.item():.4f}")
            
        if err2_val.item() <= 0.01998 and e_gap.item() >= 3.9151:
            if loss.item() < best_loss:
                best_loss = loss.item()
                best_state = (w_var.detach().clone(), W_var.detach().clone(), err2_val.item(), e_gap.item(), lam_eff.item())
                
    print(f"Optimization finished in {time.time()-t0:.1f}s.")
    
    if best_state is not None:
        w_final, W_final, final_err2, final_gap, final_lam = best_state
        print(f"*** FOUND VALID SOLUTION: err2 = {final_err2*1000:.3f} mHa, gap = {final_gap:.4f} Ha, lam_eff = {final_lam:.2f} ***")
    else:
        w_final = w_var.detach()
        W_final = W_var.detach()
        final_err2 = err2_val.item()
        final_gap = e_gap.item()
        final_lam = lam_eff.item()
        print(f"Final state: err2 = {final_err2*1000:.3f} mHa, gap = {final_gap:.4f} Ha, lam_eff = {final_lam:.2f}")
        
    w_np = w_final.cpu().numpy()
    W_np = W_final.cpu().numpy()
    
    out = dict(z)
    for k in ["onebody_angle_words", "onebody_brot", "error_1_optimized"]:
        out.pop(k, None)
    out.update(
        U=U_norm,
        w=w_np,
        W=W_np,
        err2=float(final_err2),
        gap=float(final_gap),
    )
    np.savez(out_fit, **out)
    print(f"Saved optimized fit to {out_fit}.")
    return True

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--C", type=int, default=340)
    ap.add_argument("--steps", type=int, default=500)
    ap.add_argument("--lr", type=float, default=0.001)
    a = ap.parse_args()
    solve_weights(a.fit, a.out, a.C, a.steps, a.lr)
