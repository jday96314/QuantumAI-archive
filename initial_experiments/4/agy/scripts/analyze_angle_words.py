"""Analyze the angle_words table in sample_submission/dfthcsos/submission.npz."""

import numpy as np

def main():
    sub = np.load("/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/submission.npz")
    words = sub["angle_words"]  # (3976, 75)
    print("angle_words shape:", words.shape, "dtype:", words.dtype)
    print("min word:", np.min(words), "max word:", np.max(words))
    print("Fraction of zero words:", np.mean(words == 0))
    
    # Check delta words: f_k(a) = batch[k] ^ batch[k-1]
    J = 8
    K = int(np.ceil(75 / J))
    print(f"J={J}, K={K}")
    
    total_nonzeros = 0
    total_elements = 0
    for k in range(K):
        b_curr = np.zeros((3976, J), dtype=np.int64)
        end_idx = min(75, (k + 1) * J)
        b_curr[:, :end_idx - k * J] = words[:, k * J : end_idx]
        if k == 0:
            delta = b_curr
        else:
            b_prev = np.zeros((3976, J), dtype=np.int64)
            end_prev = min(75, k * J)
            b_prev[:, :end_prev - (k - 1) * J] = words[:, (k - 1) * J : end_prev]
            delta = b_curr ^ b_prev
        nz = np.count_nonzero(delta)
        total_nonzeros += nz
        total_elements += delta.size
        print(f"Circuit {k:2d}: nonzeros = {nz}/{delta.size} ({nz/delta.size*100:.1f}%) | distinct patterns = {len(np.unique(delta, axis=0))}")
    print(f"Overall delta nonzeros: {total_nonzeros}/{total_elements} ({total_nonzeros/total_elements*100:.1f}%)")

if __name__ == "__main__":
    main()
