"""Pure linear inverse square root with F=10 (no quadratic term, no Newton refinement)."""

import sys
import math
from pathlib import Path
import numpy as np

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km

# Test floating point table generation for F=10
F = 10
h = 2.0**(-F)
num_intervals = 1 << F

print(f"Building linear coefficient table for F={F}...")
table = np.zeros((num_intervals * 2, 2), dtype=float)

for parity in [0, 1]:
    scale = 1.0 if parity == 0 else 1.0 / np.sqrt(2.0)
    for i in range(num_intervals):
        x0 = 1.0 + i * h
        x1 = x0 + h
        mid = 0.5 * (x0 + x1)
        half_w = 0.5 * h
        t0 = mid - half_w / np.sqrt(2)
        t1 = mid + half_w / np.sqrt(2)
        y0 = scale / np.sqrt(t0)
        y1 = scale / np.sqrt(t1)
        
        slope = (y1 - y0) / (t1 - t0)
        intercept = y0 - slope * t0
        
        idx = (parity << F) | i
        table[idx, 0] = intercept
        table[idx, 1] = slope

print(f"Table shape: {table.shape}")
print(f"Intercept range: [{table[:, 0].min():.4f}, {table[:, 0].max():.4f}]")
print(f"Slope range: [{table[:, 1].min():.4f}, {table[:, 1].max():.4f}]")
