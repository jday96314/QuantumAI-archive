"""Verify that our Givens tree encoder exactly matches the baseline angle_words."""

from __future__ import annotations

import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers")
import numpy as np
from dfthcsos.scoreinfra.dfthc_math import build_h1_prime, reconstruct_u_vectors
from dfthcsos.verify import native_inputs
from dfthcsos.scoreinfra.reference import check_fcidump, load_factorization, load_fixed_shifts, load_reference

def encode_givens_words(vectors: np.ndarray, tree_edges: np.ndarray, brot: int) -> np.ndarray:
    v = vectors.copy() / np.linalg.norm(vectors, axis=1, keepdims=True)
    num_vecs, num_orb = v.shape
    num_angles = num_orb - 1
    words = np.zeros((num_vecs, num_angles), dtype=np.int64)
    scale = 1 << (brot - 1)
    mod = 1 << brot

    # In reverse tree order, compute theta_t and undo the rotation
    for k in reversed(range(num_angles)):
        p, c = tree_edges[k]
        theta = np.arctan2(v[:, c], v[:, p])
        w_k = np.round((theta / np.pi) * scale).astype(np.int64) % mod
        words[:, k] = w_k
        th_q = np.pi * (w_k.astype(np.float64) / scale)
        cos_q = np.cos(th_q)
        sin_q = np.sin(th_q)
        # Inverse Givens rotation
        v[:, p] = cos_q * v[:, p] + sin_q * v[:, c]
        v[:, c] = 0.0
    return words

def main():
    sub = np.load("/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/submission.npz")
    tree_edges = np.asarray(sub["tree_edges"], dtype=np.int64)
    base_words = np.asarray(sub["angle_words"], dtype=np.int64)
    brot = 12

    # Two-body vectors from baseline U
    U = np.asarray(sub["U"], dtype=np.float64)
    R, B, N = U.shape
    u_2b = U.reshape(R * B, N)

    # Let's decode the two-body vectors from base_words[76:] and re-encode:
    u_decoded = reconstruct_u_vectors(base_words[76:] / (1 << (brot - 1)), tree_edges)
    words_2b_re = encode_givens_words(u_decoded, tree_edges, brot)
    match_2b = np.array_equal(words_2b_re, base_words[76:])
    print(f"Two-body words exact match: {match_2b}")
    if not match_2b:
        diff = np.abs(words_2b_re - base_words[76:])
        print(f"Max diff: {np.max(diff)}, num diffs: {np.sum(diff != 0)}")

    # One-body vectors from base_words[:76]
    v_decoded = reconstruct_u_vectors(base_words[:76] / (1 << (brot - 1)), tree_edges)
    words_1b_re = encode_givens_words(v_decoded, tree_edges, brot)
    match_1b = np.array_equal(words_1b_re, base_words[:76])
    print(f"One-body words exact match: {match_1b}")
    if not match_1b:
        diff = np.abs(words_1b_re - base_words[:76])
        print(f"Max diff: {np.max(diff)}, num diffs: {np.sum(diff != 0)}")

if __name__ == "__main__":
    main()
