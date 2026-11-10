"""Analyze the distribution of rank-one components across the leading eigenvectors of G."""

import numpy as np
import torch

device = torch.device("cuda")
sp = np.load("../../3/parameters/reference_spectrum.npz")
ii = sp["ii"]
jj = sp["jj"]
fac = sp["fac"]
ev = sp["eigenvalues"]
evec = sp["eigenvectors"] # shape (2926, 2926)
N = 76

# Take the top K eigenvectors of G
# Weight each eigenvector by sqrt(lambda_k)
K = 350
M_total = np.zeros((N, N))
all_evecs = []
all_weights = []

for k in range(K):
    if ev[k] <= 0:
        continue
    vk = evec[:, k]
    M_k = np.zeros((N, N))
    M_k[ii, jj] = vk / fac
    M_k[jj, ii] = vk / fac
    w, v = np.linalg.eigh(M_k) # v has shape (76, 76), w has shape (76,)
    # Energy in this eigenvector of G is ev[k] * w_j^2
    weights = np.sqrt(ev[k]) * np.abs(w)
    for j in range(N):
        if weights[j] > 1e-4:
            all_evecs.append(v[:, j])
            all_weights.append(weights[j])

all_evecs = np.array(all_evecs) # shape (M, 76)
all_weights = np.array(all_weights)
print(f"Generated {len(all_evecs)} orbital-pair eigenvectors from top {K} modes of G.")
print(f"Top 10 weights: {np.sort(all_weights)[::-1][:10].round(4)}")
print(f"Total weight sum: {np.sum(all_weights):.2f}")
