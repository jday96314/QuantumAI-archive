"""Inspect global500_lbfgs_polish.npz and positive725_c500.npz."""

from pathlib import Path
import numpy as np

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")

for name in [
    "global500_lbfgs_polish.npz",
    "global500_coordinate12.npz",
    "positive725_c500.npz",
    "residual719_lbfgs2_onebody.npz",
]:
    path = ROOT / f"initial_experiments/3/candidates/{name}"
    if not path.exists():
        continue
    z = np.load(path)
    print(f"\n--- {name} ---")
    for k in z.keys():
        val = z[k]
        if hasattr(val, "shape"):
            print(f"  {k}: shape {val.shape}, dtype {val.dtype}, min {val.min():.4g}, max {val.max():.4g}")
        else:
            print(f"  {k}: {val}")
