"""Inspect chem2b256c180_fit.npz."""

from pathlib import Path
import numpy as np

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
path = ROOT / "initial_experiments/3/candidates/chem2b256c180_fit.npz"
if path.exists():
    z = np.load(path)
    for k in z.keys():
        print(f"{k}: shape {z[k].shape}, min {z[k].min():.4g}, max {z[k].max():.4g}")
else:
    print("Does not exist")
