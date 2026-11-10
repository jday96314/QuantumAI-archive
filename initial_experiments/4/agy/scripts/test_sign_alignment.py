"""Test sign alignment of eigenvectors before Givens tree encoding."""

import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers")
import numpy as np
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
from verify_word_encoder import encode_givens_words

def main():
    sub = np.load("/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/submission.npz")
    tree_edges = np.asarray(sub["tree_edges"], dtype=np.int64)
    base_words = np.asarray(sub["angle_words"], dtype=np.int64)
    brot = 12

    # Baseline 1-body vectors
    v_base = reconstruct_u_vectors(base_words[:76] / (1 << (brot - 1)), tree_edges)

    # Re-encode v_base directly vs with sign alignment
    words_raw = encode_givens_words(v_base, tree_edges, brot)
    v_re = reconstruct_u_vectors(words_raw / (1 << (brot - 1)), tree_edges)
    err_raw = np.linalg.norm(v_re - v_base)
    print(f"v_base round-trip error: {err_raw:.6e}")

    # Now take arbitrary eigenvectors and test with/without sign alignment
    rng = np.random.default_rng(42)
    H = rng.standard_normal((76, 76))
    H = H + H.T
    _, evecs = np.linalg.eigh(H) # shape (76, 76)
    
    # 1. Raw evecs.T
    w_raw = encode_givens_words(evecs.T, tree_edges, brot)
    v_rec_raw = reconstruct_u_vectors(w_raw / (1 << (brot - 1)), tree_edges)
    # Check overlaps
    overlaps_raw = np.abs(np.diag(v_rec_raw @ evecs))
    print(f"Raw evecs mean overlap: {np.mean(overlaps_raw):.6f}, min: {np.min(overlaps_raw):.6f}")

    # 2. Check root orbital sign after tree traversal
    # In tree_edges, what is the root?
    # Orbital 0 is the root!
    # Let's align each eigenvector so its root value after reverse traversal is positive!
    v_test = evecs.T.copy()
    for k in reversed(range(75)):
        p, c = tree_edges[k]
        theta = np.arctan2(v_test[:, c], v_test[:, p])
        th_q = np.pi * (np.round((theta / np.pi) * 2048) / 2048)
        cos_q = np.cos(th_q)
        sin_q = np.sin(th_q)
        v_test[:, p] = cos_q * v_test[:, p] + sin_q * v_test[:, c]
        v_test[:, c] = 0.0
    root_signs = np.sign(v_test[:, 0])
    print("Root signs distribution:", np.unique(root_signs, return_counts=True))
    
    # Flip eigenvectors that have negative root sign!
    evecs_aligned = (evecs.T * root_signs[:, None]).copy()
    w_aligned = encode_givens_words(evecs_aligned, tree_edges, brot)
    v_rec_aligned = reconstruct_u_vectors(w_aligned / (1 << (brot - 1)), tree_edges)
    overlaps_aligned = np.abs(np.diag(v_rec_aligned @ evecs))
    diff_aligned = np.linalg.norm(v_rec_aligned - evecs_aligned, axis=1)
    print(f"Aligned evecs mean overlap: {np.mean(overlaps_aligned):.8f}, min: {np.min(overlaps_aligned):.8f}")
    print(f"Max vector L2 error with sign alignment: {np.max(diff_aligned):.6e}")

if __name__ == "__main__":
    main()
