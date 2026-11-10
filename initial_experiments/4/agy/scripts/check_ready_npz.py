#!/usr/bin/env python3
import numpy as np
from pathlib import Path

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
for name in ["residual719_q8_c350_onebody_ready.npz", "grow725_q8_c350_onebody.npz"]:
    p = ROOT / "initial_experiments/3/candidates" / name
    if p.exists():
        z = np.load(p)
        print(name, "keys:", z.files)
        words = z.get("angle_words", None)
        print(f"  words shape: {getattr(words, 'shape', None)}, brot: {z.get('brot')}")
