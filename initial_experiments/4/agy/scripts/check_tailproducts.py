#!/usr/bin/env python3
import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km
from invsqrt.verify_submission_invsqrt import verify, ProblemSpecInvSqrt

cand = ROOT / "initial_experiments/3/candidates/tailproducts_l17_q97.kmb"
print("Checking tailproducts_l17_q97.kmb...")
c = km.Circuit.from_file(cand, format="kmb")
print(f"qubits: {c.num_qubits}, toffolis: {c.operation_statistics().num_toffoli_like}")
res = verify(c, ProblemSpecInvSqrt.for_n_p(10), max_workers=8)
print("Result:", res)
