#!/usr/bin/env python3
import numpy as np
from pathlib import Path

p = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/candidates/chem35_constrained.npz")
z = np.load(p)
print("lambda_eff:", z.get("lambda_eff"))
print("err2:", z.get("err2"))
