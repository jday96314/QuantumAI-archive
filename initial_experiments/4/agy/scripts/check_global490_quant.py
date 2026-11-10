"""Check if global490_anneal_fit is quantized."""

from pathlib import Path
import numpy as np

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
path = ROOT / "initial_experiments/3/candidates/global490_anneal_fit.npz"
z = np.load(path)
print("Keys:", list(z.keys()))
if "brot" in z:
    print("brot:", z["brot"])
if "angle_words" in z:
    print("angle_words:", z["angle_words"].shape, z["angle_words"].max())
