#!/usr/bin/env python3
import numpy as np
from pathlib import Path

p = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/candidates/lbfgs500_q10c375_submission/submission.npz")
z = np.load(p)
print("Keys:", z.files)
for k in z.files:
    v = z[k]
    print(f"  {k:15s}: shape {getattr(v, 'shape', None)}, dtype {getattr(v, 'dtype', type(v))}")
