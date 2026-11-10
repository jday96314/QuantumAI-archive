"""Compare X_p.T @ X_p with reconstruct_h2(X_rho_c)."""

import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers")
import numpy as np
from dfthcsos.verify import load_inputs
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors, split_u_vectors, reconstruct_x, reconstruct_h2

def main():
    inputs = load_inputs("/tmp/c24_submission.npz")
    sub = np.load("/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/submission.npz")
    tree_edges = np.asarray(sub["tree_edges"], dtype=np.int64)
    base_words = np.asarray(sub["angle_words"], dtype=np.int64)
    brot = 12

    u_all = reconstruct_u_vectors(base_words / (1 << (brot - 1)), tree_edges)
    _, u_b_rho = split_u_vectors(u_all, inputs.w)

    X_rho_c = reconstruct_x(u_b_rho, inputs.w)
    X_rho_c = np.tril(X_rho_c) + np.swapaxes(np.tril(X_rho_c, k=-1), -1, -2)

    h2 = reconstruct_h2(X_rho_c)
    tei = inputs.ref["tei"]
    S = inputs.S
    diagonal = np.arange(76)

    h2_delta = h2 - tei
    h2_delta[:, :, diagonal, diagonal] -= 0.5 * S[:, :, None]
    h2_delta[diagonal, diagonal, :, :] -= 0.5 * S[None, :, :]
    err2 = float(np.linalg.norm(h2_delta.ravel()))
    print(f"Direct h2 err2: {err2:.6f}")

    # Now let's check baseline w on the same u_b_rho
    w_base = np.asarray(sub["W"], dtype=np.float64)
    X_base = reconstruct_x(u_b_rho, w_base)
    X_base = np.tril(X_base) + np.swapaxes(np.tril(X_base, k=-1), -1, -2)
    h2_base = reconstruct_h2(X_base)
    h2_delta_base = h2_base - tei
    h2_delta_base[:, :, diagonal, diagonal] -= 0.5 * S[:, :, None]
    h2_delta_base[diagonal, diagonal, :, :] -= 0.5 * S[None, :, :]
    err2_base = float(np.linalg.norm(h2_delta_base.ravel()))
    print(f"Baseline w err2: {err2_base:.6f}")

if __name__ == "__main__":
    main()
