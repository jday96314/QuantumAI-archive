"""Inspect all DQI stream results in initial_experiments/3/results/."""

from pathlib import Path
import json

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
results = ROOT / "initial_experiments/3/results"

base = 24439463583672.93

print(f"{'Name':<45} | {'Qubits':>6} | {'Toffoli':>8} | {'Rep':>6} | {'Cost':>18} | {'Ratio':>7}")
print("-" * 100)

for f in sorted(results.glob("stream*verify.json")):
    try:
        d = json.loads(f.read_text())
        q = d.get("num_qubits", 0)
        t = d.get("toffoli_count", 0)
        r = d.get("repetition_factor", 1.0)
        cost = d.get("score", float("nan"))
        ratio = base / cost if cost else 0
        print(f"{f.stem:<45} | {q:6d} | {t:8.0f} | {r:6.3f} | {cost:18.2f} | {ratio:6.2f}x")
    except Exception as e:
        pass
