"""Check projection floor of global500_tree_q8_0."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
import numpy as np
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

z = np.load(ROOT / "initial_experiments/3/parameters/global500_tree_q8_0.npz")
words = z["angle_words"]
tree = z["tree_edges"]
brot = int(z["brot"])

U_q = reconstruct_u_vectors(words / (1 << (brot - 1)), tree)
print(f"Reconstructed U_q shape: {U_q.shape}")

sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
ii, jj, fac = sp["ii"], sp["jj"], sp["fac"]
ev, V = sp["eigenvalues"], sp["eigenvectors"]
G = (V * np.maximum(ev, 0)) @ V.T
total_norm_sq = float(np.square(np.maximum(ev, 0)).sum())

P = U_q[:, ii] * U_q[:, jj] * fac
S = P @ P.T
se, sv = np.linalg.eigh(S)
pos = se > 1e-10
sv = sv[:, pos]
se = se[pos]
iroot = (sv / np.sqrt(se)) @ sv.T
Q = iroot @ P
H = Q @ G @ Q.T
inside_sq = float(np.square(H).sum())
outside = max(0.0, total_norm_sq - inside_sq)
floor = np.sqrt(outside)
print(f"global500_tree_q8_0 -> Quantized 8-bit projection floor = {floor*1000:.3f} mHa")
