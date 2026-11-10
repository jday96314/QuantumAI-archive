#!/usr/bin/env python3
import sys
from pathlib import Path

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path.insert(0, str(ROOT / "initial_experiments/3/vendor/uncompute"))
sys.path.insert(0, str(ROOT / "initial_experiments/2"))
import importlib.util

SOURCE = ROOT / "initial_experiments/3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/circuit.py"
spec = importlib.util.spec_from_file_location("circuit", SOURCE)
mod = importlib.util.module_from_spec(spec)
mod.ZMASKI = [0] * 70
mod.assert_metrics = lambda b, p: b._stream_totals()[0]
spec.loader.exec_module(mod)

b, toff = mod.build()
print(f"Discovered build: Toffoli = {toff}, peak_qubits = {b.peak_qubits}, next_qubit = {b.next_qubit}")
print(f"Free qubits at end: {len(b.free_qubits)}")

# Let's see the wire log or allocations in b
print("Allocated qubits total:", b.next_qubit)
