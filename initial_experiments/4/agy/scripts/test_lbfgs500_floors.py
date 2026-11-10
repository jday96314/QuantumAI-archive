import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
import numpy as np
import torch
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
from compress_chem import encode

sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
ii, jj, fac = sp["ii"], sp["jj"], sp["fac"]
ev = sp["eigenvalues"]
G = (sp["eigenvectors"] * np.maximum(ev, 0)) @ sp["eigenvectors"].T
total_norm_sq = float(np.sum(np.maximum(ev, 0)**2))

z = np.load("../../3/candidates/lbfgs500_q10c375_submission/submission.npz")
U = z["U"][0]
tree = z["tree_edges"]

# Let's check floor at brot=8, 9, 10
for brot in [8, 9, 10]:
    words = encode(U, tree, brot)
    U_q = reconstruct_u_vectors(words / (1 << (brot - 1)), tree)
    u_norm = U_q / np.linalg.norm(U_q, axis=1, keepdims=True)
    P = u_norm[:, ii] * u_norm[:, jj] * fac
    S = P @ P.T + 1e-6 * np.eye(len(U_q))
    se, sv = np.linalg.eigh(S)
    pos = se > 1e-10
    iroot = (sv[:, pos] / torch.from_numpy(np.sqrt(se[pos]))).numpy() @ sv[:, pos].T
    Q = iroot @ P
    H = Q @ G @ Q.T
    E = np.sum(H**2)
    floor = np.sqrt(max(0.0, total_norm_sq - E))
    print(f"lbfgs500 at brot={brot:2d}: floor = {floor*1000:7.4f} mHa")
