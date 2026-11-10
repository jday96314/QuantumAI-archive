"""Test post-quantization GPU optimization of float64 (w, W) at brot = 7, 8, 9, 10, 11."""

from __future__ import annotations

import time
import numpy as np
import torch
from dfthcsos.scoreinfra.dfthc_math import evaluate_dfthc_angle_words, reconstruct_u_vectors
from dfthcsos.scoreinfra.reference import UPPER_ENERGY


def encode_givens_words(vectors: np.ndarray, tree_edges: np.ndarray, brot: int) -> np.ndarray:
    """Encode unit vectors into integer Givens words in [0, 2**brot) on tree_edges."""
    v = vectors.copy() / np.linalg.norm(vectors, axis=1, keepdims=True)
    # Align sign so root 0 ends up positive after reverse rotations
    num_vecs, num_orb = v.shape
    num_angles = num_orb - 1
    words = np.zeros((num_vecs, num_angles), dtype=np.int64)
    scale = 1 << (brot - 1)
    mod = 1 << brot

    # First pass to determine root sign
    v_tmp = v.copy()
    for k in reversed(range(num_angles)):
        p, c = tree_edges[k]
        v_tmp[:, p] = np.hypot(v_tmp[:, p], v_tmp[:, c])
    # Wait: in forward pass, v starts with v[:, 0] = +1.0!
    # At each step (p, c), if we use theta = arctan2(v[:, c], v[:, p]), then before that step
    # v[:, p] = +hypot(v[:, p], v[:, c]) >= 0 (except possibly at k=0 where sign of v[:, p] matters!).
    # Wait: if k=0 is the first edge from root 0, then before k=0, v[:, 0] = +1.0 > 0, so even at k=0,
    # using arctan2(v[:, c], v[:, p]) works for ALL signs without any sign flip!
    for k in reversed(range(num_angles)):
        p, c = tree_edges[k]
        theta = np.arctan2(v[:, c], v[:, p])
        w_k = np.round((theta / np.pi) * scale).astype(np.int64) % mod
        words[:, k] = w_k
        # Apply exact inverse of the QUANTIZED rotation so parent absorbs child rounding!
        th_q = np.pi * (w_k.astype(np.float64) / scale)
        cos_q = np.cos(th_q)
        sin_q = np.sin(th_q)
        # Forward was: p_new = cos*p_old - sin*c_old, c_new = sin*p_old + cos*c_old
        # Inverse is:  p_old = cos*p_new + sin*c_new
        v[:, p] = cos_q * v[:, p] + sin_q * v[:, c]
        v[:, c] = 0.0
    return words


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
    R, B, C = w_base.shape
    N = 76

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

    for brot in [7, 8, 9, 10, 11]:
        t0 = time.time()
        # 1. Quantize U_norm at brot bits and decode exact quantized two-body vectors U_q
        words_u = encode_givens_words(U_norm.reshape(R * B, N), tree_edges, brot)
        U_q = reconstruct_u_vectors(words_u / (1 << (brot - 1)), tree_edges).reshape(R, B, N)
        U_q_t = torch.tensor(U_q, dtype=torch.float64, device=device)
        ub_q_t = (U_q_t[:, :, rows_t] * U_q_t[:, :, cols_t]) * weights_t[None, None, :]  # (65, 60, 2926)

        # Check error_2 before optimizing w, then optimize w on GPU!
        w_var = torch.tensor(w_base, dtype=torch.float64, device=device, requires_grad=True)
        W_var = torch.tensor(W_base, dtype=torch.float64, device=device, requires_grad=True)

        with torch.no_grad():
            X_p0 = torch.einsum("rbc,rbp->rcp", w_var, ub_q_t).reshape(R * C, 2926)
            err2_init = torch.norm(X_p0.T @ X_p0 - M_t).item()

        opt_w = torch.optim.LBFGS([w_var], lr=1.0, max_iter=80, history_size=30, line_search_fn="strong_wolfe")

        def closure_w() -> torch.Tensor:
            opt_w.zero_grad()
            X_p = torch.einsum("rbc,rbp->rcp", w_var, ub_q_t).reshape(R * C, 2926)
            diff = X_p.T @ X_p - M_t
            loss2 = torch.sum(diff * diff)
            # Also lightly regularize L1 norm of w so Lambda stays low
            w_B = W_var.detach() + torch.sum(w_var, dim=1)
            lam_h2 = 0.25 * torch.sum((torch.sum(torch.abs(w_var), dim=1) + torch.abs(w_B)) ** 2)
            loss = loss2 + 1e-7 * lam_h2
            loss.backward()
            return loss

        opt_w.step(closure_w)
        with torch.no_grad():
            X_p = torch.einsum("rbc,rbp->rcp", w_var, ub_q_t).reshape(R, C, 2926)
            err2_opt = torch.norm(X_p.reshape(R * C, 2926).T @ X_p.reshape(R * C, 2926) - M_t).item()

            # Build full symmetric X_rc (R, C, N, N)
            X_rc = torch.zeros((R, C, N, N), dtype=torch.float64, device=device)
            X_rc[:, :, rows_t, cols_t] = X_p / weights_t[None, None, :]
            X_rc = X_rc + X_rc.transpose(-1, -2) - torch.diagonal(X_rc, dim1=-2, dim2=-1).diag_embed()

            # Compute initial h1_prime and quantize its 76 eigenvectors to brot bits
            h1_p = h1_fixed_t - torch.einsum("rc,rcpq->pq", W_var, X_rc)
            evals0, evecs0 = torch.linalg.eigh(h1_p)
            words_1b = encode_givens_words(evecs0.T.cpu().numpy(), tree_edges, brot)
            V_q = reconstruct_u_vectors(words_1b / (1 << (brot - 1)), tree_edges)  # (N, N)
            V_q_t = torch.tensor(V_q, dtype=torch.float64, device=device)
            err1_init = torch.norm(torch.einsum("j,jp,jq->pq", evals0, V_q_t, V_q_t) - h1_p).item()

        # Now optimize W_var so h1_prime(W) aligns with the quantized eigenvectors V_q_t!
        # Wait: for fixed V_q_t, the best eigenvalues are lambda_j = diag(V_q_t @ h1_p @ V_q_t.T) (or eigh(h1_p))!
        opt_W = torch.optim.LBFGS([W_var], lr=1.0, max_iter=100, history_size=30, line_search_fn="strong_wolfe")

        def closure_W() -> torch.Tensor:
            opt_W.zero_grad()
            h1_cur = h1_fixed_t - torch.einsum("rc,rcpq->pq", W_var, X_rc)
            ev_cur = torch.linalg.eigvalsh(h1_cur)
            h1_recon = torch.einsum("j,jp,jq->pq", ev_cur, V_q_t, V_q_t)
            err1_sq = torch.sum((h1_recon - h1_cur) ** 2)
            # Keep E_gap > 3.92
            lam_h1 = torch.sum(torch.abs(ev_cur))
            e_sos = ecore - electrons * beta_1 - (lam_h1 - torch.sum(ev_cur)) - 0.5 * torch.sum(W_var ** 2)
            e_gap = energy_upper - e_sos
            gap_pen = torch.relu(3.925 - e_gap) ** 2 + torch.relu(e_gap - 3.96) ** 2
            loss = err1_sq + 10.0 * gap_pen
            loss.backward()
            return loss

        opt_W.step(closure_W)
        with torch.no_grad():
            h1_cur = h1_fixed_t - torch.einsum("rc,rcpq->pq", W_var, X_rc)
            ev_cur = torch.linalg.eigvalsh(h1_cur)
            err1_opt = torch.norm(torch.einsum("j,jp,jq->pq", ev_cur, V_q_t, V_q_t) - h1_cur).item()

        # Verify with exact evaluate_dfthc_angle_words
        all_words = np.vstack([words_1b, words_u])
        ref_dict = {"oei": cache["oei"], "tei": tei, "electrons": electrons, "ecore": ecore}
        metrics = evaluate_dfthc_angle_words(
            all_words,
            w_var.detach().cpu().numpy(),
            W_var.detach().cpu().numpy(),
            S,
            beta_1,
            ref_dict,
            brot=brot,
            tree_edges=tree_edges,
            energy_upper=energy_upper,
        )
        print(
            f"brot={brot:2d}: err2 {err2_init:.5f}->{metrics['error_2']:.5f}, "
            f"err1 {err1_init:.5f}->{metrics['error_1']:.5f}, "
            f"Lam={metrics['Lambda']:.3f}, E_gap={metrics['E_gap']:.5f}, lam_eff={metrics['lambda_eff']:.4f} "
            f"({time.time()-t0:.1f}s)",
            flush=True,
        )


if __name__ == "__main__":
    main()
