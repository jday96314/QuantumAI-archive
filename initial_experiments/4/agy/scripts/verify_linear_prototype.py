"""Verify linear_invsqrt_f10.kmb on 10,000 shots."""
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
print("Loaded circuit. Running verification on 10,000 shots...")
res = verify(c, ProblemSpecInvSqrt.for_n_p(10), max_workers=8, shots_per_x=1)
print(res)
