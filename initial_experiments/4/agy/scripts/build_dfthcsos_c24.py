"""Build dfthcsos submission with C=24 and verify with from_host/verifiers/dfthcsos/verify.py."""

from __future__ import annotations

import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers")
import json
import shutil
import time
from pathlib import Path
import numpy as np
import torch
from dfthcsos.verify import verify

def optimize_and_save_c24():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    base_sub = np.load(
        "/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/submission.npz",
        allow_pickle=False,
    )
    cache = np.load("/tmp/femoco_ref_cache.npz", allow_pickle=False)
    U_raw = np.asarray(base_sub["U"], dtype=np.float64)
    U_norm = U_raw / np.linalg.norm(U_raw, axis=2, keepdims=True)
    w_base = np.asarray(base_sub["W"], dtype=np.float64)  # (65, 60, 60)
    W_base = np.asarray(base_sub["H2_I"], dtype=np.float64) - np.sum(w_base, axis=1)
    tree_edges = np.asarray(base_sub["tree_edges"], dtype=np.int64)
    N = 76
    R, B, C_orig = w_base.shape

    S = np.asarray(cache["S"], dtype=np.float64)
    beta_1 = float(cache["beta_1"])
    h1_thc = np.asarray(cache["h1_thc"], dtype=np.float64)
    tei = np.asarray(cache["tei"], dtype=np.float64)
    ecore = float(cache["ecore"])
    electrons = int(cache["electrons"])
    energy_upper = -2452.6079299093525 # UPPER_ENERGY

    rows, cols = np.tril_indices(N)
    weights_tril = np.where(rows == cols, 1.0, np.sqrt(2.0))

    h2_target = tei.copy()
    diag = np.arange(N)
    h2_target[:, :, diag, diag] += 0.5 * S[:, :, None]
    h2_target[diag, diag, :, :] += 0.5 * S[None, :, :]
    M_pack = h2_target[rows[:, None], cols[:, None], rows[None, :], cols[None, :]]
    M_pack *= weights_tril[:, None] * weights_tril[None, :]

    M_t = torch.tensor(M_pack, dtype=torch.float64, device=device)
    weights_t = torch.tensor(weights_tril, dtype=torch.float64, device=device)
    rows_t = torch.tensor(rows, dtype=torch.long, device=device)
    cols_t = torch.tensor(cols, dtype=torch.long, device=device)
    h1_fixed = h1_thc - 0.5 * electrons * S + beta_1 * np.eye(N)
    h1_fixed_t = torch.tensor(h1_fixed, dtype=torch.float64, device=device)

    U_t = torch.tensor(U_norm, dtype=torch.float64, device=device)
    U_b_tril = (U_t[:, :, rows_t] * U_t[:, :, cols_t]) * weights_t[None, None, :]

    w_t = torch.tensor(w_base, dtype=torch.float64, device=device)
    W_t = torch.tensor(W_base, dtype=torch.float64, device=device)
    w_rot = torch.zeros_like(w_t)
    W_rot = torch.zeros_like(W_t)
    for r in range(R):
        xr = torch.einsum("bc,bp->cp", w_t[r], U_b_tril[r])
        Q_r, _, _ = torch.linalg.svd(xr, full_matrices=False)
        w_rot[r] = w_t[r] @ Q_r
        W_rot[r] = W_t[r] @ Q_r

    C_new = 24
    w_var = w_rot[:, :, :C_new].clone().detach().requires_grad_(True)
    W_var = W_rot[:, :C_new].clone().detach().requires_grad_(True)

    opt = torch.optim.LBFGS([w_var, W_var], lr=1.0, max_iter=200, history_size=40, line_search_fn="strong_wolfe")

    def closure() -> torch.Tensor:
        opt.zero_grad()
        X_p = torch.einsum("rbc,rbp->rcp", w_var, U_b_tril).reshape(R * C_new, 2926)
        diff = X_p.T @ X_p - M_t
        loss2 = torch.sum(diff * diff)
        X_rc = torch.zeros((R, C_new, N, N), dtype=torch.float64, device=device)
        X_rc[:, :, rows_t, cols_t] = X_p.reshape(R, C_new, 2926) / weights_t[None, None, :]
        X_rc = X_rc + X_rc.transpose(-1, -2) - torch.diagonal(X_rc, dim1=-2, dim2=-1).diag_embed()
        h1_p = h1_fixed_t - torch.einsum("rc,rcpq->pq", W_var, X_rc)
        ev = torch.linalg.eigvalsh(h1_p)
        lam_h1 = torch.sum(torch.abs(ev))
        e_sos = ecore - electrons * beta_1 - (lam_h1 - torch.sum(ev)) - 0.5 * torch.sum(W_var ** 2)
        e_gap = energy_upper - e_sos
        gap_penalty = torch.relu(3.93 - e_gap) ** 2 + torch.relu(e_gap - 4.1) ** 2
        loss = loss2 + 1e-4 * gap_penalty
        loss.backward()
        return loss

    print("Optimizing w and W for C=24 on GPU...")
    t0 = time.time()
    opt.step(closure)
    print(f"Optimization finished in {time.time()-t0:.1f}s.")

    w_opt = w_var.detach().cpu().numpy()
    W_opt = W_var.detach().cpu().numpy()
    H2_I_opt = W_opt + np.sum(w_opt, axis=1)

    out_dir = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/dfthcsos")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_npz = out_dir / "submission.npz"

    np.savez(
        out_npz,
        U=base_sub["U"],
        W=w_opt,
        H2_I=H2_I_opt,
        H1_I=base_sub["H1_I"],
        Hsym_lower=base_sub["Hsym_lower"],
        tree_edges=tree_edges,
    )
    print(f"Saved {out_npz} with C=24 (W shape: {w_opt.shape}, H2_I shape: {H2_I_opt.shape})")

    # Copy angle circuits if not present
    for k in range(10):
        src = f"/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/angle_k_{k}.kmb"
        dst = out_dir / f"angle_k_{k}.kmb"
        if not dst.exists():
            shutil.copy(src, dst)

def main():
    optimize_and_save_c24()
    out_dir = "/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/dfthcsos"
    print("\nRunning official verify() on new dfthcsos/...")
    res = verify(out_dir, format="kmb")
    print("Verification result:")
    for k, v in res.items():
        if k not in ("resources", "metrics"):
            print(f"  {k}: {v}")

    with open("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/logs/verify_dfthcsos.json", "w") as f:
        json.dump(res, f, indent=2)
    print("Saved logs/verify_dfthcsos.json")

if __name__ == "__main__":
    main()
