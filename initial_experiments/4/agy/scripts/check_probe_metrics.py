#!/usr/bin/env python3
import numpy as np
from pathlib import Path

p = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/candidates/global500_onebody_q8_probe.npz")
z = np.load(p)
print("err2:", float(z["err2"])*1000, "mHa")
print("err1:", float(z["error_1_optimized"])*1000, "mHa")
print("brot:", z["brot"])
