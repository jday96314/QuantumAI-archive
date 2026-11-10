#!/usr/bin/env python3
import numpy as np
from pathlib import Path

p = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/candidates/global500_onebody_q8_probe.npz")
z = np.load(p)
print("Keys:", z.files)
for k in z.files:
    v = z[k]
    print(f"  {k:20s}: shape {getattr(v, 'shape', None)}, dtype {getattr(v, 'dtype', type(v))}")
