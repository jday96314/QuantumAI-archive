#!/usr/bin/env python3
import numpy as np
from pathlib import Path

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
for name in ["chem35_constrained.npz", "chem40_constrained.npz", "chem12b128c90_fit.npz"]:
    p = ROOT / "initial_experiments/3/candidates" / name
    if p.exists():
        z = np.load(p)
        print(name, "keys:", z.files)
        U = z["U"]
        w = z["w"]
        err2 = z.get("err2", None)
        print(f"  U: {U.shape}, w: {w.shape}, err2: {err2}")
