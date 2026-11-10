"""Analyze the matrix structure of the leading eigenvectors of G."""

import numpy as np

sp = np.load("../../3/parameters/reference_spectrum.npz")
ii, jj, fac = sp["ii"], sp["jj"], sp["fac"]
ev = sp["eigenvalues"]
evec = sp["eigenvectors"] # shape (2926, 2926)

print(f"Top 10 eigenvalues of G:")
for k in range(10):
    print(f"  k={k}: lambda = {ev[k]:.4f}")

# Look at the first 5 eigenvectors as 76x76 matrices
N = 76
for k in [0, 1, 2, 5, 10, 50, 100]:
    vk = evec[:, k]
    M_k = np.zeros((N, N))
    M_k[ii, jj] = vk / fac
    M_k[jj, ii] = vk / fac # symmetric
    w, _ = np.linalg.eigh(M_k)
    w_abs = np.sort(np.abs(w))[::-1]
    top3 = w_abs[:5]
    print(f"Eigenvector k={k:3d} (lambda={ev[k]:8.3f}): top 5 singular values of 76x76 matrix M_k: {top3.round(3)}, trace = {np.trace(M_k):.3f}")
