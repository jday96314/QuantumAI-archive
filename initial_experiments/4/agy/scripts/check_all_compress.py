#!/usr/bin/env python3
import numpy as np
from pathlib import Path

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
for c in [347, 348, 349, 350]:
    p = ROOT / f"initial_experiments/3/parameters/residual719_compress_c{c}_init.npz"
    if p.exists():
        z = np.load(p)
        print(f"c={c}: err2={float(z['err2'])*1000:.3f} mHa, brot={z.get('brot')}")
