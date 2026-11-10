"""Test exact brot=11 with fixed Givens encoding on baseline factorization."""

from __future__ import annotations

import numpy as np
import kickmix as km
from dfthcsos.scoreinfra.dfthc_math import evaluate_dfthc_angle_words
from dfthcsos.scoreinfra.reference import UPPER_ENERGY


def encode_givens(vecs: np.ndarray, edges: np.ndarray, brot: int) -> np.ndarray:
    v = vecs.copy()
    A, N = v.shape
    words = np.zeros((A, N - 1), dtype=np.int64)
    scale = 1 << (brot - 1)
    mask = (1 << brot) - 1
    for k in reversed(range(N - 1)):
        p, c = edges[k]
        vp = v[:, p]
        vc = v[:, c]
        hyp = np.hypot(vp, vc)
        theta = np.arctan2(vc, vp)
        theta = np.where(theta < 0, theta + 2 * np.pi, theta)
        w = np.round(theta * (scale / np.pi)).astype(np.int64) & mask
        words[:, k] = w
        v[:, p] = hyp
        v[:, c] = 0.0
    return words


def main() -> None:
    sub = np.load("/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/submission.npz")
    cache = np.load("/tmp/femoco_ref_cache.npz")
    edges = sub["tree_edges"]

    # 1-body vectors from h1_prime of baseline
    w = np.asarray(sub["W"], dtype=np.float64)
    W = np.asarray(sub["H2_I"], dtype=np.float64) - np.sum(w, axis=1)
    S = np.asarray(cache["S"], dtype=np.float64)
    beta_1 = float(cache["beta_1"])

    U_npz = np.asarray(sub["U"], dtype=np.float64)
    U_norm = U_npz / np.linalg.norm(U_npz, axis=2, keepdims=True)
    u_2body = U_norm.reshape(-1, 76)

    # Let's get the 1-body vectors from baseline submission angle_words
    base_words = sub["angle_words"]
    # Reconstruct vectors from base_words
    from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
    vecs_all = reconstruct_u_vectors(base_words / (1 << 11), edges)
    u_1body = vecs_all[:76]

    v_all = np.vstack([u_1body, u_2body])
    print(f"Total vectors: {v_all.shape}")

    for b in [12, 11]:
        w_enc = encode_givens(v_all, edges, b)
        metrics = evaluate_dfthc_angle_words(
            w_enc, w, W, S, beta_1, cache, brot=b, tree_edges=edges, energy_upper=UPPER_ENERGY
        )
        print(f"brot={b}: error_1={metrics['error_1']*1000:.3f} mHa, error_2={metrics['error_2']*1000:.3f} mHa, lambda_eff={metrics['lambda_eff']:.4f}")


if __name__ == "__main__":
    main()
