"""Check err2 and C of union700_q8_c340.npz."""

import numpy as np
from pathlib import Path

p = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/candidates/union700_q8_c340.npz")
z = np.load(p)
print("err2:", float(z['err2']) * 1000, "mHa")
if "w" in z:
    print("w shape:", z["w"].shape)
if "K" in z:
    print("K shape:", z["K"].shape)
