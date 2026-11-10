"""Inspect union700_q8_c340.npz."""

from pathlib import Path
import numpy as np

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
path = ROOT / "initial_experiments/3/candidates/union700_q8_c340.npz"
z = np.load(path)
print("Keys:", list(z.keys()))
for k in z.keys():
    val = z[k]
    if hasattr(val, "shape"):
        print(f"{k}: shape {val.shape}, dtype {val.dtype}, min {val.min():.4g}, max {val.max():.4g}")
    else:
        print(f"{k}: {val}")
