#!/usr/bin/env python3
import numpy as np
from pathlib import Path

p = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/residual719_compress_c350_init.npz")
z = np.load(p)
print("Keys:", z.files)
for k in z.files:
    v = z[k]
    print(f"  {k:15s}: shape {getattr(v, 'shape', None)}, dtype {getattr(v, 'dtype', type(v))}")
