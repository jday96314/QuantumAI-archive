"""Check stats of all DQI candidate circuits in experiment 3 candidates."""

from pathlib import Path
import sys
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km

base = 24439463583672.93

candidates_dir = ROOT / "initial_experiments/3/candidates"
for p in sorted(candidates_dir.glob("dqi*.kmb")):
    try:
        c = km.Circuit.from_file(p, format="kmb")
        q = c.num_qubits
        t = c.operation_statistics().num_toffoli_like
        est_cost = 1.22 * (q ** 2) * t
        ratio = base / est_cost
        print(f"{p.name:25s} | Q={q:4d} | T={t:10,d} | Est Cost: {est_cost:16,.2f} ({ratio:6.2f}x)")
    except Exception as e:
        print(f"{p.name:25s} | Error: {e}")
