#!/usr/bin/env python3
import numpy as np
from pathlib import Path

p = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/residual719_compress_c350_init.npz")
z = np.load(p)
print("err2:", float(z["err2"])*1000, "mHa")
print("brot:", z["brot"])
