"""Run greedy dictionary selection for B=250 and B=300 from existing vectors."""

import sys
import json
import time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
ev, V = sp["eigenvalues"], sp["eigenvectors"]
ii, jj = sp["ii"], sp["jj"]
fac = sp["fac"]

# Pool of vectors from chem35 and residual719
z1 = np.load(ROOT / "initial_experiments/3/candidates/chem35c45_rot.npz")
u1 = reconstruct_u_vectors(z1["angle_words"][76:] / 2048, z1["tree_edges"])

z2 = np.load(ROOT / "initial_experiments/3/candidates/residual719_lbfgs2_onebody.npz")
u2 = reconstruct_u_vectors(z2["angle_words"][76:] / 128, z2["tree_edges"])

U_pool = np.vstack([u1, u2])
print(f"Total candidate vectors in pool: {len(U_pool)}")

P = U_pool[:, ii] * U_pool[:, jj] * fac
F = (V[:, :600] * np.sqrt(np.maximum(ev[:600], 0))).T

gram = (U_pool @ U_pool.T) ** 2
cross = P @ F.T

for B in [350, 300, 250]:
    start = time.time()
    g = gram.copy()
    cr = cross.copy()
    chosen = []
    for k in range(B):
        diag = np.maximum(g.diagonal(), 1e-20)
        score = (cr * cr).sum(1) / diag
        score[chosen] = -1
        j = int(score.argmax())
        chosen.append(j)
        v = g[:, j] / np.sqrt(diag[j])
        h = cr[j] / np.sqrt(diag[j])
        g -= v[:, None] * v[None, :]
        cr -= v[:, None] * h[None, :]

    u_sel = U_pool[chosen]
    P_sel = P[chosen]
    s = P_sel @ P_sel.T
    se, sv = eigh(s)
    se = np.maximum(se, 1e-10)
    iroot = (sv / np.sqrt(se)) @ sv.T
    Q = iroot @ P_sel
    G = (V * ev) @ V.T
    H = Q @ G @ Q.T
    e, v = eigh(H)
    e = np.maximum(e, 0)
    order = np.argsort(e)[::-1]
    e = e[order]
    v = v[:, order]
    w = (iroot @ v) * np.sqrt(e)
    min_err2 = float(np.sqrt(max(0, np.square(ev).sum() - np.square(e).sum())))
    print(f"B={B:3d} | Selected in {time.time()-start:.1f}s | min_error_2: {min_err2*1000:.3f} mHa")
    
    # Save candidate
    out_path = f"/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/runs/isdf_b{B}.npz"
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    np.savez(out_path, U=u_sel[None, :, :], w=w[None, :, :], W=np.zeros((1, B)), chosen=np.array(chosen), min_err2=min_err2)
    print(f"Saved {out_path}")
