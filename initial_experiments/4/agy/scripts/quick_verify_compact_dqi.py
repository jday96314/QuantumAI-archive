"""Quick verification of dqi_compact_psi.kmb on 256 shots."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km
from dqi.verify_submission_dqi import verify, ProblemSpecDQI

c_path = ROOT / "initial_experiments/4/agy/runs/dqi_compact_psi.kmb"
print(f"Loading {c_path}...")
c = km.Circuit.from_file(c_path, format="kmb")
print(f"Loaded: {c.num_qubits} qubits, {c.operation_statistics().num_toffoli_like:,} Toffolis.")

print("Running verification on 256 shots with 4 workers...")
spec = ProblemSpecDQI.benchmark_instance(num_shots=256)
res = verify(c, spec, max_workers=4)

print("\nVerification result:")
for k, v in res.items():
    if k != "shot_operation_counts":
        print(f"  {k}: {v}")

base = 24439463583672.93
score = res.get("score")
if score:
    print(f"\nCost: {score:,.2f} | Speedup: {base / score:.4f}x")
