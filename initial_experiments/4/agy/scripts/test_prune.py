"""Prune global500 down to smaller basis sizes and observe error_2."""

import subprocess
import sys
from pathlib import Path

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")

fit_path = ROOT / "initial_experiments/3/candidates/global500_lbfgs_polish.npz"
out_dir = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/runs")
out_dir.mkdir(parents=True, exist_ok=True)

out_file = out_dir / "pruned_global350.npz"

cmd = [
    sys.executable,
    str(ROOT / "initial_experiments/3/prune_global.py"),
    "--fit", str(fit_path),
    "--out", str(out_file),
    "--basis", "350",
]

print("Running command:", " ".join(cmd))
subprocess.run(cmd, check=True)
