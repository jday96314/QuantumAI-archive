#!/usr/bin/env python3
import json
from pathlib import Path

hook_path = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/stream10_offset3_select_noscale_compact_hook.json")
with open(hook_path) as f:
    hook = json.load(f)

print("maxq:", hook.get("maxq"))
print("free wires count:", len(hook.get("free", [])))
print("free wires:", hook.get("free", [])[:20])
print("sigma count:", len(hook.get("sigma", [])))
print("tail count:", len(hook.get("tail", [])))

# Let's see the range of wire IDs used in stream10
all_sigma = [w for reg in hook.get("sigma", []) for w in reg]
all_tail = [w for reg in hook.get("tail", []) for w in reg]
print("Min sigma wire:", min(all_sigma), "Max sigma wire:", max(all_sigma))
print("Min tail wire:", min(all_tail), "Max tail wire:", max(all_tail))
