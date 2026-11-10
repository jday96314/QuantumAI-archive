"""Inspect positive700_c335_tree.npz and its quantized counterpart."""

import numpy as np
from pathlib import Path

p = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/positive700_c335_tree.npz")
z = np.load(p)
print("Keys:", list(z.keys()))
for k in z.keys():
    arr = z[k]
    if isinstance(arr, np.ndarray):
        print(f"  {k}: shape {arr.shape}, dtype {arr.dtype}")
    else:
        print(f"  {k}: {arr}")

pq = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/positive700_c335_tree_quantized.npz")
zq = np.load(pq)
print("\nQuantized Keys:", list(zq.keys()))
for k in zq.keys():
    arr = zq[k]
    if isinstance(arr, np.ndarray):
        print(f"  {k}: shape {arr.shape}, dtype {arr.dtype}")
    else:
        print(f"  {k}: {arr}")
