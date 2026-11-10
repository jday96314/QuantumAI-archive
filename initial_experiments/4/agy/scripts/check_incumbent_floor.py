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

sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
ii, jj, fac = sp["ii"], sp["jj"], sp["fac"]
ev = sp["eigenvalues"]
G = (sp["eigenvectors"] * np.maximum(ev, 0)) @ sp["eigenvectors"].T
total_norm_sq = float(np.sum(np.maximum(ev, 0)**2))

sub = np.load("dfthcsos/submission.npz")
U = sub["U"][0]
tree = sub["tree_edges"]
words = sub["angle_words"].reshape(-1, 75)
U_q = reconstruct_u_vectors(words / 128.0, tree)
u_norm = U_q / np.linalg.norm(U_q, axis=1, keepdims=True)
P = u_norm[:, ii] * u_norm[:, jj] * fac
S = P @ P.T + 1e-6 * np.eye(len(U_q))
se, sv = np.linalg.eigh(S)
pos = se > 1e-10
iroot = (sv[:, pos] / np.sqrt(se[pos])) @ sv[:, pos].T
Q = iroot @ P
H = Q @ G @ Q.T
E = np.sum(H**2)
floor = np.sqrt(max(0.0, total_norm_sq - E))
print(f"Incumbent 719 atoms quantized floor: {floor*1000:.4f} mHa")
