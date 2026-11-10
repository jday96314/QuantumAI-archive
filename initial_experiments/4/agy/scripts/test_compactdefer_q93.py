#!/usr/bin/env python3
"""Test candidate compactdefer_k3_q93.kmb on host invsqrt verifier."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
import kickmix as km
from invsqrt.verify_submission_invsqrt import verify, ProblemSpecInvSqrt

path = ROOT / "initial_experiments/3/candidates/compactdefer_k3_q93.kmb"
print(f"Loading {path.name}...")
c = km.Circuit.from_file(path, format="kmb")
print(f"qubits: {c.num_qubits}, toffolis: {c.operation_statistics().num_toffoli_like}")

# First quick test on smaller shot batch to check validity
spec_small = ProblemSpecInvSqrt(n_p=10, num_shots=100000, batch_size=256)
res = verify(c, spec_small, max_workers=8)
print("Small verification result:", res)
