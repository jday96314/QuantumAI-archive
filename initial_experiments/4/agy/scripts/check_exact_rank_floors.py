"""Compute exact subspace projection error floor vs B for R=1, 2, 3, 4."""

import numpy as np

ROOT = Path = "/media/shared/ExpansionDrive2/QuantumAI"
sp = np.load(ROOT + "/initial_experiments/3/parameters/reference_spectrum.npz")
ev, V = sp["eigenvalues"], sp["eigenvectors"]
ii, jj, fac = sp["ii"], sp["jj"], sp["fac"]
G = (V * np.maximum(ev, 0)) @ V.T
total_norm_sq = float(np.square(np.maximum(ev, 0)).sum())

# SVD of G: G has shape (2926, 2926).
# The eigenvalues of G give the exact optimal low-rank approximation in Frobenius norm!
print("Spectrum of G:")
sorted_ev = np.sort(np.maximum(ev, 0))[::-1]
for rank in [100, 150, 200, 250, 300, 350, 400, 500, 719]:
    # Best rank-k approximation error:
    # ||G - G_k||_F = sqrt(sum_{i > k} lambda_i^2)
    err = np.sqrt(np.sum(sorted_ev[rank:]**2))
    print(f"Rank {rank:4d} optimal error floor: {err*1000:7.3f} mHa")
