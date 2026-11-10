"""Inspect hook JSON of stream10."""

import json
from pathlib import Path

hook_path = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/stream10_offset3_select_noscale_compact_hook.json")
d = json.loads(hook_path.read_text())
print("Keys in hook:", list(d.keys()))
for k in d.keys():
    if k not in ("phases", "sigma", "free", "omega", "eo", "ab"):
        print(f"  {k}: {d[k]}")
print(f"free len: {len(d.get('free', []))}")
print(f"sigma len: {len(d.get('sigma', []))}")
