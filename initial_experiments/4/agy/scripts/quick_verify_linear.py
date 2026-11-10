"""Quick test of linear_invsqrt_f10.kmb on small n_p."""
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km
from invsqrt.verify_submission_invsqrt import verify, ProblemSpecInvSqrt

c_path = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/runs/linear_invsqrt_f10.kmb")
c = km.Circuit.from_file(c_path, format="kmb")
print("Loaded circuit. Running verification on n_p=6 (quick sanity check)...")
spec = ProblemSpecInvSqrt.for_n_p(10)
print(f"Running full verification (200,540,928 shots) on {spec.max_x} values...")
res = verify(c, spec, max_workers=8)
print("Result:", res)
