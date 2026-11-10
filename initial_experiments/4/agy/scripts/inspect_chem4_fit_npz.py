#!/usr/bin/env python3
import numpy as np
from pathlib import Path

cand = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/candidates/chem4b256c128_fit.npz")
z = np.load(cand)
print("Keys in chem4b256c128_fit.npz:", z.files)
for k in z.files:
    v = z[k]
    print(f"  {k:15s}: shape {getattr(v, 'shape', None)}, dtype {getattr(v, 'dtype', type(v))}")
