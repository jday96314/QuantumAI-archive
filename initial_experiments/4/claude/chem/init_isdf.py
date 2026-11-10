"""Greedy ISDF-style block dictionaries from a candidate vector pool.

The leading R*C pair-space eigen-directions of the target are dealt to R blocks
(interleaved by eigenvalue); each block greedily selects B candidate vectors that
best capture its own target directions (pivoted Gram elimination).
"""
import argparse, sys, numpy as np
sys.path.insert(0, '/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers')
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

p = argparse.ArgumentParser()
p.add_argument('--R', type=int, required=True); p.add_argument('--B', type=int, required=True)
p.add_argument('--C', type=int, required=True); p.add_argument('--out', required=True)
p.add_argument('--directions', type=int, default=600)
a = p.parse_args()
E3 = '/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/'
sp = np.load(E3 + 'parameters/reference_spectrum.npz')
ev, V, ii, jj, fac = sp['eigenvalues'], sp['eigenvectors'], sp['ii'], sp['jj'], sp['fac']
order = np.argsort(ev)[::-1]
pool = []
z = np.load(E3 + 'candidates/chem35c45_rot.npz')
pool.append(reconstruct_u_vectors(z['angle_words'][76:] / 2048, z['tree_edges']))
pool.append(np.load(E3 + 'candidates/mixing_mu1e7.npz')['U'][0])
n = 76
eig = []
for c in order[:a.directions]:
    m = np.zeros((n, n)); m[ii, jj] = V[:, c] / fac; m = m + m.T - np.diag(np.diag(m))
    mu, vec = np.linalg.eigh(m)
    top = np.argsort(np.abs(mu))[::-1][:6]
    eig.append(vec[:, top].T)
pool.append(np.concatenate(eig))
U = np.concatenate(pool); U /= np.linalg.norm(U, axis=1, keepdims=True)
print('pool', U.shape, flush=True)
P = U[:, ii] * U[:, jj] * fac
dirs = order[:a.directions]
groups = [dirs[r::a.R] for r in range(a.R)]
used = np.zeros(len(U), bool)
blocks = []
gram_all = (U @ U.T) ** 2
for r, g in enumerate(groups):
    F = (V[:, g] * np.sqrt(np.maximum(ev[g], 0))).T
    gmat = gram_all.copy(); cr = P @ F.T; chosen = []
    for k in range(a.B):
        diag = np.maximum(gmat.diagonal(), 1e-20); score = (cr * cr).sum(1) / diag
        score[used] = -1; score[chosen] = -1
        j = int(score.argmax()); chosen.append(j)
        v = gmat[:, j] / np.sqrt(diag[j]); h = cr[j] / np.sqrt(diag[j])
        gmat -= v[:, None] * v[None, :]; cr -= v[:, None] * h[None, :]
    used[chosen] = True; blocks.append(U[chosen])
    print('block', r, 'captured', float(1 - (cr * cr).sum() / (F * F).sum()), flush=True)
np.savez(a.out, U=np.stack(blocks), w=np.zeros((a.R, a.B, a.C)))
