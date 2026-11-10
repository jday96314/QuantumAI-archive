"""Analyze qubit allocations across phases in stream10."""
import json
from pathlib import Path
import numpy as np

hook = json.loads(Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/stream10_offset3_select_noscale_compact_hook.json").read_text())
print("Total ops in tape:", hook["ops"])
print("maxq:", hook["maxq"])
print("peak:", hook["peak"])
print("outcomes:", hook["outcomes"])

phases = hook.get("phases", [])
print(f"Total marked phases: {len(phases)}")
for idx, name in phases[:30]:
    print(f"  Op {idx:8d}: {name}")
