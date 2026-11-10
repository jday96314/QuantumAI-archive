"""Fast SVD analysis of w[rho, b, c] and X^(rho, c) in baseline dfthcsos using cached integrals."""

from __future__ import annotations

import numpy as np
from dfthcsos.verify import heuristic_costs
from dfthcsos.scoreinfra.dfthc_math import lambda_h2, symmetric_from_lower
from dfthcsos.scoreinfra.reference import UPPER_ENERGY


def main() -> None:
    sub = np.load(
        "/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/submission.npz",
        allow_pickle=False,
    )
    cache = np.load("/tmp/femoco_ref_cache.npz", allow_pickle=False)
    U_raw = np.asarray(sub["U"], dtype=np.float64)
    U_norm = U_raw / np.linalg.norm(U_raw, axis=2, keepdims=True)
    w = np.asarray(sub["W"], dtype=np.float64)  # (65, 60, 60)
    W = np.asarray(sub["H2_I"], dtype=np.float64) - np.sum(w, axis=1)  # (65, 60)
    R, B, C = w.shape
    N = 76

    S = np.asarray(cache["S"], dtype=np.float64)
    beta_1 = float(cache["beta_1"])
    h1_thc = np.asarray(cache["h1_thc"], dtype=np.float64)
    tei = np.asarray(cache["tei"], dtype=np.float64)
    ecore = float(cache["ecore"])
    electrons = int(cache["electrons"])
    energy_upper = float(UPPER_ENERGY)

    # Packed 2926 lower triangle with Frobenius weights (1 on diag, sqrt(2) off-diag)
    rows, cols = np.tril_indices(N)
    weights_tril = np.where(rows == cols, 1.0, np.sqrt(2.0))

    # Target M_pack (2926, 2926)
    h2_target = tei.copy()
    diag = np.arange(N)
    h2_target[:, :, diag, diag] += 0.5 * S[:, :, None]
    h2_target[diag, diag, :, :] += 0.5 * S[None, :, :]
    M_pack = h2_target[rows[:, None], cols[:, None], rows[None, :], cols[None, :]]
    M_pack *= weights_tril[:, None] * weights_tril[None, :]

    # For each rho, compute U_b_tril (B, 2926) and X_r (C, 2926)
    U_b_all = np.zeros((R, B, len(rows)), dtype=np.float64)
    sv_c_all = np.zeros((R, C))
    Q_all = np.zeros((R, C, C))
    for r in range(R):
        u_r = U_norm[r]
        U_b_tril = (u_r[:, rows] * u_r[:, cols]) * weights_tril[None, :]
        U_b_all[r] = U_b_tril
        X_r = w[r].T @ U_b_tril
        Q_r, s_c, _ = np.linalg.svd(X_r, full_matrices=False)
        sv_c_all[r] = s_c
        Q_all[r] = Q_r

    print("Top 10 mean singular values across C (per rho):", np.round(np.mean(sv_c_all[:, :10], axis=0), 6))
    print("Top 10 max singular values across C (per rho):", np.round(np.max(sv_c_all[:, :10], axis=0), 6))

    for c_new in [20, 22, 24, 25, 26, 28, 30, 35, 40, 50]:
        w_new = np.einsum("rbc,rck->rbk", w, Q_all[:, :, :c_new])
        W_new = np.einsum("rc,rck->rk", W, Q_all[:, :, :c_new])
        # X_pack: (R * c_new, 2926)
        X_pack = np.einsum("rbk,rbp->rkp", w_new, U_b_all).reshape(R * c_new, -1)
        err2 = float(np.linalg.norm(X_pack.T @ X_pack - M_pack, ord="fro"))

        # Reconstruct sum_{r,c} W_rc X_rc
        WX_pack = np.einsum("rk,rkp->p", W_new, X_pack.reshape(R, c_new, -1)) / weights_tril
        WX_mat = np.zeros((N, N), dtype=np.float64)
        WX_mat[rows, cols] = WX_pack
        WX_mat = WX_mat + WX_mat.T - np.diag(np.diag(WX_mat))

        h1_prime = h1_thc - 0.5 * electrons * S + beta_1 * np.eye(N) - WX_mat
        evals = np.linalg.eigvalsh(h1_prime)
        lam_h1 = float(np.sum(np.abs(evals)))
        Lam = lam_h1 + lambda_h2(w_new, W_new)
        E_SOS = ecore - electrons * beta_1 - (lam_h1 - float(np.sum(evals))) - 0.5 * float(np.sum(W_new**2))
        E_gap = energy_upper - E_SOS
        lam_eff = np.sqrt(E_gap * (2.0 * Lam - E_gap)) if E_gap > 0 else float("nan")
        h_cost = heuristic_costs(N, R, B, c_new, 11)
        print(
            f"C_new={c_new:2d}: err2={err2:.6f}, Lam={Lam:.4f}, E_gap={E_gap:.6f}, lam_eff={lam_eff:.4f} | "
            f"H_T={h_cost['toffoli']}, H_Q={h_cost['qubits']} (temp={h_cost['temporary_qubits']})",
            flush=True,
        )


if __name__ == "__main__":
    main()
