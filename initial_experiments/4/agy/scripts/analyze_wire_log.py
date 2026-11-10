"""Analyze wire_log in Builder to trace live wires and peak occupancy."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "initial_experiments/3/vendor/uncompute"),
]
import importlib.util

spec = importlib.util.spec_from_file_location(
    "discovered_circuit",
    ROOT / "initial_experiments/3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/circuit.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

mod.ZMASKI = [0] * 70
b, t = mod.build()

print("Analyzing wire_log events:", len(b.wire_log))
print("Sample event 0:", b.wire_log[0])
print("Sample event 1:", b.wire_log[1])

live_by_set = {}
peak_live = 0
peak_breakdown = {}

current_live = 0
for item in b.wire_log:
    ev = item[0]
    if ev == "claim":
        name = item[1][0]
        ws = item[1][1]
        cat = name.split("/")[0]
        live_by_set[cat] = live_by_set.get(cat, 0) + len(ws)
        current_live += len(ws)
        if current_live > peak_live:
            peak_live = current_live
            peak_breakdown = dict(live_by_set)
    elif ev == "release":
        name = item[1][0]
        ws = item[1][1]
        cat = name.split("/")[0]
        live_by_set[cat] = live_by_set.get(cat, 0) - len(ws)
        current_live -= len(ws)

print(f"Recorded Peak live wires: {peak_live}")
print("Breakdown at peak:")
for k, v in sorted(peak_breakdown.items(), key=lambda x: x[1], reverse=True):
    print(f"  {k:20s}: {v:5d} qubits ({v/12:.1f} words)")
