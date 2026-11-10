"""Inspect stages and resource breakdown of circuit.py in vendor/uncompute."""

from pathlib import Path
import sys

CIRCUIT_PY = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/circuit.py")
content = CIRCUIT_PY.read_text().splitlines()

for i, line in enumerate(content):
    if "def build" in line or "STAGE" in line or "Stage" in line or "stage" in line and "def" in line:
        print(f"Line {i+1}: {line}")
