"""Inspect the exact recovered angle words and vectors of sample_submission/dfthcsos."""

from __future__ import annotations

import numpy as np
from dfthcsos.scoreinfra.circuits import load_circuits, recover_words_from_circuits
from dfthcsos.scoreinfra.dfthc_math import evaluate_dfthc_vectors, reconstruct_u_vectors
from dfthcsos.scoreinfra.reference import UPPER_ENERGY


def main() -> None:
    sub_dir = "/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos"
    sub = np.load(f"{sub_dir}/submission.npz", allow_pickle=False)
    cache = np.load("/tmp/femoco_ref_cache.npz", allow_pickle=False)

    circuits = load_circuits(sub_dir)
    words, layout, toff = recover_words_from_circuits(
        circuits, num_vectors=3976, num_angles=75, shots_per_address=1, seed=0
    )
    print(f"Baseline circuits: brot={layout.brot}, J={layout.batch_size}, toff={toff}")

    tree_edges = np.asarray(sub["tree_edges"], dtype=np.int64)
    u_from_circuits = reconstruct_u_vectors(words / (1 << (layout.brot - 1)), tree_edges)

    U_npz = np.asarray(sub["U"], dtype=np.float64)
    U_npz_norm = U_npz / np.linalg.norm(U_npz, axis=2, keepdims=True)
    u_2b_circuits = u_from_circuits[76:].reshape(65, 60, 76)

    # Check overlap between u_2b_circuits and U_npz_norm!
    dots = np.abs(np.sum(u_2b_circuits * U_npz_norm, axis=2))
    print("Min/Mean/Max |dot| between circuit 2-body vectors and U_npz_norm:", np.min(dots), np.mean(dots), np.max(dots))

    # Wait: what is U in submission.npz?! Could U in submission.npz be angle_words or unnormalized or transposed or something else?!
    w = np.asarray(sub["W"], dtype=np.float64)
    W = np.asarray(sub["H2_I"], dtype=np.float64) - np.sum(w, axis=1)
    S = np.asarray(cache["S"], dtype=np.float64)
    beta_1 = float(cache["beta_1"])
    ref_dict = {"oei": cache["oei"], "tei": cache["tei"], "electrons": int(cache["electrons"]), "ecore": float(cache["ecore"])}

    m_circ = evaluate_dfthc_vectors(u_from_circuits, w, W, S, beta_1, ref_dict, energy_upper=float(UPPER_ENERGY))
    print("Metrics from baseline circuits:", m_circ)

    # Check how words[76:] in baseline circuits relate to U_npz_norm!
    def encode_givens(vecs: np.ndarray, edges: np.ndarray, brot: int) -> np.ndarray:
        v = vecs / np.linalg.norm(vecs, axis=1, keepdims=True)
        out = np.zeros((len(v), len(edges)), dtype=np.int64)
        scale = 1 << (brot - 1)
        mod = 1 << brot
        for k in reversed(range(len(edges))):
            p, c = edges[k]
            th = np.arctan2(v[:, c], v[:, p])
            out[:, k] = np.round((th / np.pi) * scale).astype(np.int64) % mod
            v[:, p] = np.hypot(v[:, p], v[:, c])
            v[:, c] = 0.0
        return out

    words_our12 = encode_givens(U_npz_norm.reshape(3900, 76), tree_edges, 12)
    u_our12 = reconstruct_u_vectors(words_our12 / 2048.0, tree_edges)
    dots_our = np.abs(np.sum(u_our12 * U_npz_norm.reshape(3900, 76), axis=1))
    print("Min/Mean/Max |dot| of our encode_givens(12) vs U_npz_norm:", np.min(dots_our), np.mean(dots_our), np.max(dots_our))
    diff_words = np.abs(words[76:] - words_our12)
    print("Max word difference between baseline circuits and our encode_givens(12):", np.max(diff_words))


if __name__ == "__main__":
    main()
