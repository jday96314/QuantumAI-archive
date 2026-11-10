"""Test angle quantization on chem4b256c128."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
import numpy as np
from compress_chem import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors, evaluate_dfthc_vectors

z = np.load(ROOT / "initial_experiments/3/candidates/chem4b256c128_float_rot_off.npz")
U = z["U"]
w = z["w"]
W = z["W"]
R, B, N = U.shape

tree = np.load(ROOT / "initial_experiments/3/candidates/chem35c45_rot.npz")["tree_edges"]
ref = np.load(ROOT / "initial_experiments/3/parameters/chem_reference.npz")
S = ref["S"]
beta_1 = float(ref["beta_1"])

print(f"chem4b256c128: R={R}, B={B}, N={N}")

for brot in [8, 9, 10, 11, 12]:
    # Encode U
    words = encode(U.reshape(-1, N), tree, brot)
    U_q = reconstruct_u_vectors(words / (1 << (brot - 1)), tree)
    
    # Evaluate X
    U_reshaped = U_q.reshape(R, B, N)
    u_norm = U_reshaped / np.linalg.norm(U_reshaped, axis=2, keepdims=True)
    X = np.einsum("tjp,tjq,tjl->tlpq", u_norm, u_norm, w, optimize=True)
    X = np.tril(X) + np.swapaxes(np.tril(X, k=-1), -1, -2)
    
    # Evaluate error_2
    sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
    ii, jj, fac = sp["ii"], sp["jj"], sp["fac"]
    ev, V = sp["eigenvalues"], sp["eigenvectors"]
    G = (V * ev) @ V.T
    
    # Flatten X to shape (R*C, 2926)
    X_flat = X.reshape(R * 128, N, N)
    P = (X_flat[:, ii, jj] * fac).reshape(R * 128, -1)
    diff = P.T @ P - G
    err2 = float(np.sqrt(np.sum(diff**2)))
    
    print(f"brot={brot:2d} -> error_2: {err2*1000:7.3f} mHa")
