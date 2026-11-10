#!/usr/bin/env python3
import json
from pathlib import Path

liveness_path = Path("../../3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/liveness.json")
with open(liveness_path) as f:
    data = json.load(f)

phase_names = [p.get("name") for p in data.get("phases", [])]
idx_euclid_R = phase_names.index("euclid/R")
phase_to_idx = {name: i for i, name in enumerate(phase_names)}

euclid_R_qubits = set()
sets_in_euclid_R = {}

for item in data.get("intervals", []):
    bp = item.get("born_phase")
    dp = item.get("died_phase")
    b_idx = phase_to_idx.get(bp, 0)
    d_idx = phase_to_idx.get(dp, len(phase_names))
    if b_idx <= idx_euclid_R <= d_idx:
        q = item.get("q")
        euclid_R_qubits.add(q)
        for s in item.get("sets", []):
            sets_in_euclid_R[s] = sets_in_euclid_R.get(s, 0) + 1

print(f"Total live qubits during euclid/R: {len(euclid_R_qubits)}")
print("Sets active during euclid/R:", sorted(sets_in_euclid_R.items(), key=lambda x: -x[1]))
