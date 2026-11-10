"""Test global SVD and GPU optimization of (U, w, W) for compressed (R, B, C) in dfthcsos."""

from __future__ import annotations

import time
import numpy as np
import torch


def main() -> None:
    device = torch.device("cuda")
    sub = np.load(
        "/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/submission.npz",
        allow_pickle=False,
    )
    cache = np.load("/tmp/femoco_ref_cache.npz", allow_pickle=False)
    U_raw = np.asarray(sub["U"], dtype=np.float64)
    U_norm = U_raw / np.linalg.norm(U_raw, axis=2, keepdims=True)
    w_base = np.asarray(sub["W"], dtype=np.float64)  # (65, 60, 60)
    W_base = np.asarray(sub["H2_I"], dtype=np.float64) - np.sum(w_base, axis=1)
    N = 76

    S = np.asarray(cache["S"], dtype=np.float64)
    beta_1 = float(cache["beta_1"])
    h1_thc = np.asarray(cache["h1_thc"], dtype=np.float64)
    tei = np.asarray(cache["tei"], dtype=np.float64)
    ecore = float(cache["ecore"])
    electrons = int(cache["electrons"])
    energy_upper = -2452.6079299093525

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

    # Check global singular values of X_base (3900 x 2926)
    U_t = torch.tensor(U_norm, dtype=torch.float64, device=device)
    w_t = torch.tensor(w_base, dtype=torch.float64, device=device)
    U_b_tril = (U_t[:, :, rows_t] * U_t[:, :, cols_t]) * weights_t[None, None, :]  # (65, 60, 2926)
    X_base = torch.einsum("rbc,rbp->rcp", w_t, U_b_tril).reshape(-1, 2926)  # (3900, 2926)
    s_global = torch.linalg.svdvals(X_base)
    for k in [200, 250, 300, 325, 350, 400, 500, 600, 800, 1000]:
        tail = torch.sqrt(torch.sum(s_global[k:] ** 4)).item()
        print(f"Global X_base^T X_base rank-{k:4d} tail norm: {tail:.6f}")

    # Test GPU optimization of (U, w, W) when pruning (R, B, C) from baseline!
    # First, sort b for each rho by ||w[r, b, :]||_2, and sort rho by Frobenius norm
    b_norms = torch.norm(w_t, dim=2)  # (65, 60)
    b_idx = torch.argsort(b_norms, dim=1, descending=True)
    U_sorted = torch.gather(U_t, 1, b_idx.unsqueeze(-1).expand(-1, -1, N))
    w_sorted = torch.gather(w_t, 1, b_idx.unsqueeze(-1).expand(-1, -1, 60))

    # For each rho, rotate C by SVD of X_r
    for r in range(65):
        u_r = U_sorted[r]
        ub_r = (u_r[:, rows_t] * u_r[:, cols_t]) * weights_t[None, :]
        xr = w_sorted[r].T @ ub_r
        Q_r, _, _ = torch.linalg.svd(xr, full_matrices=False)
        w_sorted[r] = w_sorted[r] @ Q_r

    rho_norms = torch.zeros(65, dtype=torch.float64, device=device)
    for r in range(65):
        u_r = U_sorted[r]
        ub_r = (u_r[:, rows_t] * u_r[:, cols_t]) * weights_t[None, :]
        xr = w_sorted[r].T @ ub_r
        rho_norms[r] = torch.norm(xr)
    rho_idx = torch.argsort(rho_norms, descending=True)
    U_sorted = U_sorted[rho_idx]
    w_sorted = w_sorted[rho_idx]

    # Now test joint optimization of (U, w) to minimize ||X^T X - M_pack||_F for compressed (R, B, C)
    for R_try, B_try, C_try in [
        (43, 16, 11),  # R*C*B = 7568 <= 8192 -> k2=4! A = 76 + 688 = 764 (10-bit selection!)
        (32, 24, 10),  # R*C*B = 7680 <= 8192 -> k2=4! A = 76 + 768 = 844 (10-bit selection!)
        (55, 24, 6),   # R*C*B = 7920 <= 8192 -> k2=4! A = 76 + 1320 = 1396 (11-bit selection!)
        (60, 32, 16),  # R*C*B = 30720 (k2=5), A = 1996 (11-bit selection)
        (65, 32, 20),  # R*C*B = 41600 (k2=5), A = 2156
    ]:
        t0 = time.time()
        U_var = U_sorted[:R_try, :B_try, :].clone().detach().requires_grad_(True)
        w_var = w_sorted[:R_try, :B_try, :C_try].clone().detach().requires_grad_(True)

        opt = torch.optim.LBFGS(
            [U_var, w_var],
            lr=1.0,
            max_iter=150,
            history_size=40,
            line_search_fn="strong_wolfe",
        )

        def closure() -> torch.Tensor:
            opt.zero_grad()
            u_n = U_var / torch.norm(U_var, dim=2, keepdim=True)
            ub = (u_n[:, :, rows_t] * u_n[:, :, cols_t]) * weights_t[None, None, :]
            X_p = torch.einsum("rbc,rbp->rcp", w_var, ub).reshape(R_try * C_try, 2926)
            diff = X_p.T @ X_p - M_t
            loss = torch.sum(diff * diff)
            loss.backward()
            return loss

        opt.step(closure)
        err2 = torch.sqrt(closure()).item()
        print(
            f"(R={R_try:2d}, B={B_try:2d}, C={C_try:2d}) [RCB={R_try*B_try*C_try:5d}, A={76+R_try*B_try:4d}]: "
            f"after 150 L-BFGS steps -> err2 = {err2:.6f} ({time.time() - t0:.1f}s)",
            flush=True,
        )


if __name__ == "__main__":
    main()
