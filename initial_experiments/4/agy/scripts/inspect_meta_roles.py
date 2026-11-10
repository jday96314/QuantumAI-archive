#!/usr/bin/env python3
import json
from pathlib import Path

meta_path = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/meta.json")
with open(meta_path) as f:
    meta = json.load(f)

print("Roles in meta.json:")
for r in meta.get("roles", []):
    print(f"  {r.get('set', 'transient'):10s}: {r.get('name')}")
