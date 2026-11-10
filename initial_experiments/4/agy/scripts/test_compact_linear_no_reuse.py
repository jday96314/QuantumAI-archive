"""Test compact linear invsqrt with reuse_normalized_top=False."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys_path = ROOT / "initial_experiments/3"
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(sys_path),
    str(ROOT / "initial_experiments/1"),
]
import kickmix as km
import subprocess
import binary_centered_invsqrt as polynomial
from invsqrt.verify_submission_invsqrt import verify, ProblemSpecInvSqrt

spec = ProblemSpecInvSqrt.for_n_p(10)

for qlim in [100, 95, 90, 85, 80, 78, 75]:
    try:
        c = polynomial.build(
            10, (28, 18), bias=0, rounded=True, booth=True,
            defer_workspace=True, borrow_integration=True, carry_bits=1,
            reuse_normalized_top=False, qubit_limit=qlim
        )
        raw = Path(f"/tmp/test_linear_no_reuse_{qlim}.raw.kmx")
        compact = Path(f"/tmp/test_linear_no_reuse_{qlim}.compact.kmx")
        c.to_file(str(raw), format="kmx")
        res = subprocess.run([str(sys_path / 'compact_wires'), str(raw), str(compact), '22'], capture_output=True, text=True, check=True)
        c_comp = km.Circuit.from_file(str(compact))
        q = c_comp.num_qubits
        t = c_comp.operation_statistics().num_toffoli_like
        est_cost = q * (t + 250)
        print(f"qubit_limit={qlim:2d} -> alloc={c.num_qubits}, compacted={q}, Toffoli={t}, est_cost={est_cost:,}")
    except Exception as e:
        print(f"qubit_limit={qlim:2d} -> FAILED: {e}")
