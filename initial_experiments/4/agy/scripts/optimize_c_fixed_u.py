"""Optimize w and W with fixed baseline U for compressed C on GPU."""

from __future__ import annotations

import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers")
import time
import numpy as np
import torch
from dfthcsos.scoreinfra.reference import UPPER_ENERGY
from dfthcsos.scoreinfra.costs_prepare import SpinFreeData_V1

def heuristic_costs(num_orbitals, rank_r, bases_b, ensemble_c, brot):
    model = SpinFreeData_V1(R=int(rank_r), B=int(bases_b), C=int(ensemble_c), N=int(num_orbitals), chi=10, beta=int(brot))
    outer_qubits = model.outer_state_prep_qubit_cost()
    inner_qubits = model.inner_state_prep_qubit_cost()
    outer_persistent = int(outer_qubits["persistent"])
    outer_temporary = int(outer_qubits["temporary"])
    inner_persistent = int(inner_qubits["persistent"])
    inner_temporary = int(inner_qubits["temporary"])
    toff = int(model.outer_state_prep() + model.inner_state_prep() + model.givens_rotations())
    q = outer_persistent + inner_persistent + max(outer_temporary, inner_temporary)
    return {"toffoli": toff, "qubits": q}


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
    tree_edges = np.asarray(sub["tree_edges"], dtype=np.int64)
    N = 76
    R, B, C_orig = w_base.shape

    S = np.asarray(cache["S"], dtype=np.float64)
    beta_1 = float(cache["beta_1"])
    h1_thc = np.asarray(cache["h1_thc"], dtype=np.float64)
    tei = np.asarray(cache["tei"], dtype=np.float64)
    ecore = float(cache["ecore"])
    electrons = int(cache["electrons"])
    energy_upper = float(UPPER_ENERGY)

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

    # Compute exact U_b_tril
    U_t = torch.tensor(U_norm, dtype=torch.float64, device=device)
    U_b_tril = (U_t[:, :, rows_t] * U_t[:, :, cols_t]) * weights_t[None, None, :]  # (65, 60, 2926)

    # Rotate baseline w into principal components per rho
    w_t = torch.tensor(w_base, dtype=torch.float64, device=device)
    W_t = torch.tensor(W_base, dtype=torch.float64, device=device)
    w_rot = torch.zeros_like(w_t)
    W_rot = torch.zeros_like(W_t)
    for r in range(R):
        xr = torch.einsum("bc,bp->cp", w_t[r], U_b_tril[r])
        Q_r, _, _ = torch.linalg.svd(xr, full_matrices=False)
        w_rot[r] = w_t[r] @ Q_r
        W_rot[r] = W_t[r] @ Q_r

    # Test compression for C in [16, 18, 20, 22, 24, 26]
    for C_new in [16, 18, 20, 22, 24, 26]:
        w_var = w_rot[:, :, :C_new].clone().detach().requires_grad_(True)
        W_var = W_rot[:, :C_new].clone().detach().requires_grad_(True)

        opt = torch.optim.LBFGS([w_var, W_var], lr=1.0, max_iter=200, history_size=40, line_search_fn="strong_wolfe")

        def closure() -> torch.Tensor:
            opt.zero_grad()
            X_p = torch.einsum("rbc,rbp->rcp", w_var, U_b_tril).reshape(R * C_new, 2926)
            diff = X_p.T @ X_p - M_t
            loss2 = torch.sum(diff * diff)
            # Reconstruct h1_prime and regularize gap and err1
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

        t0 = time.time()
        opt.step(closure)
        with torch.no_grad():
            X_p = torch.einsum("rbc,rbp->rcp", w_var, U_b_tril).reshape(R * C_new, 2926)
            err2 = torch.norm(X_p.T @ X_p - M_t).item()
            X_rc = torch.zeros((R, C_new, N, N), dtype=torch.float64, device=device)
            X_rc[:, :, rows_t, cols_t] = X_p.reshape(R, C_new, 2926) / weights_t[None, None, :]
            X_rc = X_rc + X_rc.transpose(-1, -2) - torch.diagonal(X_rc, dim1=-2, dim2=-1).diag_embed()
            h1_p = h1_fixed_t - torch.einsum("rc,rcpq->pq", W_var, X_rc)
            ev = torch.linalg.eigvalsh(h1_p)
            lam_h1 = torch.sum(torch.abs(ev))
            w_B = W_var + torch.sum(w_var, dim=1)
            lam_h2 = 0.25 * torch.sum((torch.sum(torch.abs(w_var), dim=1) + torch.abs(w_B)) ** 2)
            Lam = (lam_h1 + lam_h2).item()
            e_sos = (ecore - electrons * beta_1 - (lam_h1 - torch.sum(ev)) - 0.5 * torch.sum(W_var ** 2)).item()
            e_gap = energy_upper - e_sos
            lam_eff = np.sqrt(e_gap * (2 * Lam - e_gap))
            h_cost = heuristic_costs(N, R, B, C_new, 12)
            print(
                f"C={C_new:2d}: err2={err2:.6f}, Lam={Lam:.3f}, E_gap={e_gap:.5f}, lam_eff={lam_eff:.4f} | "
                f"H_T={h_cost['toffoli']}, H_Q={h_cost['qubits']} ({time.time()-t0:.1f}s)",
                flush=True,
            )


if __name__ == "__main__":
    main()
