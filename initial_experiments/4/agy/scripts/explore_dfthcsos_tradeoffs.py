"""Explore (R, B, C) pruning/compression and cost tradeoffs for dfthcsos on GPU."""

from __future__ import annotations

import math
import numpy as np
import torch


def qr_opt(L: int, M1: int) -> tuple[int, int]:
    k = max(0.0, 0.5 * math.log2(L / M1))
    k_fl, k_ce = int(math.floor(k)), int(math.ceil(k))
    v_fl = math.ceil(L / (2**k_fl) + M1 * ((2**k_fl) - 1))
    v_ce = math.ceil(L / (2**k_ce) + M1 * ((2**k_ce) - 1))
    return (k_fl, v_fl) if v_fl <= v_ce else (k_ce, v_ce)


def qi_opt(L: int) -> tuple[int, int]:
    return qr_opt(L, 1)


def heuristic_cost_fast(R: int, B: int, C: int, brot: int = 11, N: int = 76, chi: int = 10) -> tuple[int, int, int, int]:
    X_o = N + R * C
    b1 = math.ceil(math.log2(X_o)) + chi
    outer_t = 8 * math.ceil(math.log2(X_o)) + qr_opt(X_o, b1)[1] + b1 + qi_opt(X_o)[1] + b1 + 1
    k_o = qr_opt(X_o, b1)[0]
    outer_p = (b1 + chi) + (math.ceil(math.log2(X_o)) + 2)
    outer_tmp = (math.ceil(math.log2(X_o / (2**k_o))) + (2**k_o) * b1) + max(0, math.ceil(math.log2(X_o)) - 2)

    b2 = math.ceil(math.log2(B)) + chi + 1
    k2 = qr_opt(R * C * B, b2)[0]
    for j in range(k2 + 1):
        kn = k2 - j
        if math.ceil(math.log2(R * C * math.ceil(B / (2**kn)))) + ((2**kn) - 1) * b2 < brot * (N - 1):
            k2 = kn
            break
    k2_t = k2
    if B < 2**k2_t:
        for _ in range(k2_t):
            k2_t -= 1
            if B >= 2**k2_t:
                break
    RC = R * C
    for_qroam = RC * math.ceil((B + 1) / (2**k2_t)) + (2**k2_t) * b2
    for_qroam_dag = qr_opt(RC, B + 1)[1] + (B + 1) + N + b2
    inner_t = for_qroam + for_qroam_dag + 8 * math.ceil(math.log2(B + 1))

    inner_p = (b2 + chi + math.ceil(math.log2(R)) + 2) + (math.ceil(math.log2(B + 1)) + 2)
    inner_tmp = max(0, math.ceil(math.log2(B + 1)) - 2) + math.ceil(math.log2(RC * math.ceil((B + 1) / (2**k2)))) + ((2**k2) - 1) * b2
    givens_t = 2 * (N - 1) * (2 * brot - 1)
    return outer_t + inner_t + givens_t, outer_p + inner_p, max(outer_tmp, inner_tmp), k2


def main() -> None:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
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

    # Check importance of each (rho, b) basis vector and each rho
    rows, cols = np.tril_indices(N)
    weights_tril = np.where(rows == cols, 1.0, np.sqrt(2.0))
    U_b_all = (U_norm[:, :, rows] * U_norm[:, :, cols]) * weights_tril[None, None, :]  # (65, 60, 2926)

    # Check singular values of w_base[r] across B (which basis vectors b are small?)
    b_norms = np.linalg.norm(w_base, axis=2)  # (65, 60)
    b_sorted = np.sort(b_norms, axis=1)[:, ::-1]
    print("Mean b-vector L2 norm at rank b=0, 10, 20, 30, 40, 50, 59:")
    for idx in [0, 10, 20, 30, 31, 40, 45, 50, 55, 59]:
        print(f"  b={idx:2d}: mean={np.mean(b_sorted[:, idx]):.6f}, max={np.max(b_sorted[:, idx]):.6f}")

    # Check rho importance: Frobenius norm of X_r = w_base[r].T @ U_b_all[r]
    rho_norms = np.zeros(65)
    for r in range(65):
        X_r = w_base[r].T @ U_b_all[r]
        rho_norms[r] = np.linalg.norm(X_r, ord="fro")
    rho_sorted_idx = np.argsort(rho_norms)[::-1]
    print("Min and max rho Frobenius norms:", np.min(rho_norms), np.max(rho_norms))

    # Print cost table for candidate (R, B, C, brot, J)
    print("\nCandidate (R, B, C, brot) total costs (assuming One-Way QROAM + del_lookup):")
    for R, B, C, brot in [
        (65, 60, 60, 11),
        (65, 60, 30, 11),
        (65, 60, 22, 11),
        (65, 48, 22, 11),
        (65, 32, 22, 11),
        (60, 32, 22, 11),
        (65, 16, 8, 11),
        (43, 16, 11, 11),
        (32, 16, 16, 11),
        (25, 20, 16, 11),
    ]:
        h_t, h_p, h_tmp, k2 = heuristic_cost_fast(R, B, C, brot)
        A = N + R * B
        sel_bits = math.ceil(math.log2(A))
        # With J=1 (75 circuits) or J=2 (38 circuits), block size M=4:
        # Each circuit: init_lookup(ceil(A/4), sel_bits-2) + 3*(J*brot) cswaps + del_lookup(A, sel_bits)
        for J in [1, 2]:
            K = math.ceil(75 / J)
            M = 4 if sel_bits + (M_try := 4) * J * brot + (sel_bits - 2) <= 128 else 2
            init_t = math.ceil(A / M) - 1
            cswap_t = (M - 1) * J * brot
            # del_lookup over A entries costs min_{k}(ceil(A/2^k) + 2^k - 1)
            del_t = qi_opt(1 << sel_bits)[1]
            ang_t_2pass = 2 * K * (init_t + cswap_t + del_t)
            T_tot = h_t + ang_t_2pass
            Q_tot = h_p + J * brot + h_tmp
            score_est = 49.02 * T_tot * Q_tot
            ratio = 3355599588.3155375 / score_est
            print(
                f"  R={R:2d}, B={B:2d}, C={C:2d}, brot={brot}, J={J} (k2={k2}, M={M}): "
                f"H_T={h_t:5d}, Ang_T={ang_t_2pass:5d}, T_tot={T_tot:5d}, Q_tot={Q_tot:3d} -> "
                f"score~{score_est/1e6:6.1f}M ({ratio:5.2f}x)"
            )


if __name__ == "__main__":
    main()
