"""Test C=28 using exact decoded u_b_rho to build h1_prime."""

import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers")
import numpy as np
from dfthcsos.scoreinfra.dfthc_math import build_h1_prime, evaluate_dfthc_angle_words, reconstruct_u_vectors, split_u_vectors, reconstruct_x
from dfthcsos.scoreinfra.reference import UPPER_ENERGY
from verify_word_encoder import encode_givens_words

def main():
    sub = np.load("/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/submission.npz")
    cache = np.load("/tmp/femoco_ref_cache.npz")
    tree_edges = np.asarray(sub["tree_edges"], dtype=np.int64)
    base_words = np.asarray(sub["angle_words"], dtype=np.int64)
    brot = 12
    N = 76

    w = np.asarray(sub["W"], dtype=np.float64)  # (65, 60, 60)
    W = np.asarray(sub["H2_I"], dtype=np.float64) - np.sum(w, axis=1)  # (65, 60)
    R, B, C = w.shape

    S = np.asarray(cache["S"], dtype=np.float64)
    beta_1 = float(cache["beta_1"])
    h1_thc = np.asarray(cache["h1_thc"], dtype=np.float64)
    tei = np.asarray(cache["tei"], dtype=np.float64)
    ecore = float(cache["ecore"])
    electrons = int(cache["electrons"])
    energy_upper = float(UPPER_ENERGY)

    rows, cols = np.tril_indices(N)
    weights_tril = np.where(rows == cols, 1.0, np.sqrt(2.0))

    # Exact decoded u_b_rho from baseline two-body words!
    u_all = reconstruct_u_vectors(base_words / (1 << (brot - 1)), tree_edges)
    _, u_b_rho = split_u_vectors(u_all, w)

    # Compute Q_r using the exact decoded u_b_rho!
    Q_all = np.zeros((R, C, C))
    for r in range(R):
        u_r = u_b_rho[r]
        U_b_tril = (u_r[:, rows] * u_r[:, cols]) * weights_tril[None, :]
        X_r = w[r].T @ U_b_tril
        Q_r, s_c, _ = np.linalg.svd(X_r, full_matrices=False)
        Q_all[r] = Q_r

    for c_new in [28, 30, 32, 35, 40, 50]:
        w_new = np.einsum("rbc,rck->rbk", w, Q_all[:, :, :c_new])
        W_new = np.einsum("rc,rck->rk", W, Q_all[:, :, :c_new])

        # Exact X_rho_c
        X_rho_c = reconstruct_x(u_b_rho, w_new)
        X_rho_c = np.tril(X_rho_c) + np.swapaxes(np.tril(X_rho_c, k=-1), -1, -2)

        # Exact build_h1_prime
        h1_prime = build_h1_prime(h1_thc, S, beta_1, electrons, W_new, X_rho_c)
        evals, evecs = np.linalg.eigh(h1_prime)

        # Encode exact eigenvectors into Givens words
        words_1b_new = encode_givens_words(evecs.T, tree_edges, brot)
        all_words_new = np.vstack([words_1b_new, base_words[76:]])

        ref_dict = {"oei": cache["oei"], "tei": tei, "electrons": electrons, "ecore": ecore}
        metrics = evaluate_dfthc_angle_words(
            all_words_new,
            w_new,
            W_new,
            S,
            beta_1,
            ref_dict,
            brot=brot,
            tree_edges=tree_edges,
            energy_upper=energy_upper,
        )
        print(f"c_new={c_new:2d}: error_1={metrics['error_1']*1000:.3f} mHa, error_2={metrics['error_2']*1000:.3f} mHa, "
              f"Lambda={metrics['Lambda']:.3f}, E_gap={metrics['E_gap']:.5f}, lam_eff={metrics['lambda_eff']:.4f}, "
              f"gap_in_domain={metrics['gap_in_domain']}")

if __name__ == "__main__":
    main()
