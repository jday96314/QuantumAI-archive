"""Check quantized projection floors of all candidate dictionaries with B <= 650."""

import glob, os, sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]
import numpy as np
import torch
from compress_chem import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

device = torch.device("cuda")
sp = np.load("../../3/parameters/reference_spectrum.npz")
ii = torch.tensor(sp["ii"], device=device)
jj = torch.tensor(sp["jj"], device=device)
fac = torch.tensor(sp["fac"], dtype=torch.float64, device=device)
ev = torch.tensor(sp["eigenvalues"], dtype=torch.float64, device=device)
V = torch.tensor(sp["eigenvectors"], dtype=torch.float64, device=device)
G = (V * torch.clamp(ev, min=0.0)) @ V.T
total_norm_sq = torch.sum(torch.clamp(ev, min=0.0)**2).item()

tree = np.load("../../3/candidates/chem35c45_rot.npz")["tree_edges"]

files = sorted(glob.glob("../../3/candidates/*.npz"))
print(f"Checking {len(files)} candidates...")

results = []
for f in files:
    try:
        z = np.load(f)
        if "U" not in z:
            continue
        U = z["U"]
        if U.ndim == 3:
            R, B, N = U.shape
        elif U.ndim == 2:
            B, N = U.shape
            R = 1
        else:
            continue
        if B > 650 or R > 1:
            continue
        
        # Check if angle_words present
        if "angle_words" in z and int(z.get("brot", 8)) == 8:
            words = z["angle_words"].reshape(-1, 75)
            # take last B words
            if len(words) >= B:
                words = words[-B:]
                U_q = reconstruct_u_vectors(words / 128.0, tree)
            else:
                U_q = U.reshape(B, N)
        else:
            # quantize to brot=8
            words = encode(U.reshape(B, N), tree, 8)
            U_q = reconstruct_u_vectors(words / 128.0, tree)
            
        U_t = torch.tensor(U_q / np.linalg.norm(U_q, axis=1, keepdims=True), dtype=torch.float64, device=device)
        P = U_t[:, ii] * U_t[:, jj] * fac
        S = P @ P.T
        se, sv = torch.linalg.eigh(S)
        pos = se > 1e-10
        sv = sv[:, pos]
        se = se[pos]
        iroot = (sv / torch.sqrt(se)) @ sv.T
        Q = iroot @ P
        H = Q @ G @ Q.T
        in_sq = torch.sum(H**2).item()
        floor = np.sqrt(max(0.0, total_norm_sq - in_sq))
        base_name = os.path.basename(f)
        results.append((floor * 1000, B, base_name))
    except Exception as e:
        pass

results.sort(key=lambda x: x[0])
print(f"\nTop 15 candidates by brot=8 projection floor:")
for fl, B, name in results[:15]:
    print(f"  Floor: {fl:7.3f} mHa | B = {B:3d} | {name}")
