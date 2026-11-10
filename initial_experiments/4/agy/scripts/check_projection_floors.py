"""Compute the subspace projection error floor (sqrt(outside)) for various candidate dictionaries."""

import sys
from pathlib import Path
import numpy as np

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
ev, V = sp["eigenvalues"], sp["eigenvectors"]
ii, jj, fac = sp["ii"], sp["jj"], sp["fac"]
G = (V * np.maximum(ev, 0)) @ V.T
total_norm_sq = float(np.square(np.maximum(ev, 0)).sum())

candidates = [
    ("residual719", ROOT / "initial_experiments/3/candidates/residual719_lbfgs2_onebody.npz"),
    ("global500_polish", ROOT / "initial_experiments/3/candidates/global500_lbfgs_polish.npz"),
    ("global500_coordinate12", ROOT / "initial_experiments/3/candidates/global500_coordinate12.npz"),
    ("chem35c45", ROOT / "initial_experiments/3/candidates/chem35c45_rot.npz"),
    ("chem4b256", ROOT / "initial_experiments/3/candidates/chem4b256c128_float_rot_off.npz"),
    ("chem2b256", ROOT / "initial_experiments/3/candidates/chem2b256c180_fit.npz"),
    ("positive725", ROOT / "initial_experiments/3/candidates/positive725_c500.npz"),
]

for name, path in candidates:
    z = np.load(path)
    U = z["U"]
    if U.ndim == 3 and U.shape[0] == 1:
        u = U[0]
        P = u[:, ii] * u[:, jj] * fac
    else:
        # Multi-rank or 3D
        U_flat = U.reshape(-1, 76)
        P = U_flat[:, ii] * U_flat[:, jj] * fac
    
    # Projection of G onto span of P
    S = P @ P.T
    se, sv = np.linalg.eigh(S)
    pos = se > 1e-12
    sv = sv[:, pos]
    se = se[pos]
    iroot = (sv / np.sqrt(se)) @ sv.T
    Q = iroot @ P
    H = Q @ G @ Q.T
    inside_sq = float(np.square(H).sum())
    outside = max(0.0, total_norm_sq - inside_sq)
    floor = np.sqrt(outside)
    print(f"{name:25s} | B_total={len(P):4d} | Projection floor sqrt(outside) = {floor*1000:7.3f} mHa")
