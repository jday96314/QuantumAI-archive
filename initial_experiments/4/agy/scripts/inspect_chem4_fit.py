"""Inspect chem4b256c128_fit.npz."""

import numpy as np
from pathlib import Path

p = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/candidates/chem4b256c128_fit.npz")
z = np.load(p)
print("Keys:", list(z.keys()))
for k in z.keys():
    arr = z[k]
    if isinstance(arr, np.ndarray):
        print(f"  {k}: shape {arr.shape}, dtype {arr.dtype}")
    else:
        print(f"  {k}: {arr}")
