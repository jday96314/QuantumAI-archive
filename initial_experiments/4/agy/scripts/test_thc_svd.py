import numpy as np

# Load reference spectrum
sp = np.load('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/reference_spectrum.npz')
ev = sp['eigenvalues']
V = sp['eigenvectors']
ii, jj = sp['ii'], sp['jj']
fac = sp['fac']

# Sort descending
order = np.argsort(ev)[::-1]
ev_sorted = ev[order]
V_sorted = V[:, order]

K = 300 # Top 300 components
tail_norm = np.sqrt(np.sum(np.maximum(ev_sorted[K:], 0)**2))
print(f"Top {K} SVD tail norm: {tail_norm * 1000:.3f} mHa")

# Each eigenvector v_k in R^2926 corresponds to a symmetric 76x76 matrix M_k
N = 76
matrices = np.zeros((K, N, N))
for k in range(K):
    vec = V_sorted[:, k] * np.sqrt(ev_sorted[k])
    M = np.zeros((N, N))
    M[ii, jj] = vec * fac
    M[jj, ii] = vec * fac
    matrices[k] = M

print(f"Constructed {K} matrices of shape {matrices.shape}")
