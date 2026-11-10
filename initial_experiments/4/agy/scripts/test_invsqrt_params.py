"""Test parameter variations around build_compact_invsqrt."""

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
import numpy as np
import subprocess
import binary_centered_invsqrt as poly
from invsqrt.verify_submission_invsqrt import verify, ProblemSpecInvSqrt

table = np.load(sys_path / 'candidates/margin_curve8_t5.9e-08.npy')
poly.table_for = lambda *args: table

spec = ProblemSpecInvSqrt.for_n_p(10)

for oeb in [10, 9, 8]:
    for plim in [17, 16, 15]:
        for qlim in [94, 90, 85]:
            try:
                c = poly.build(
                    7, (27, 18, 8), lookup_frac=8, accum_bits=[18, 9, None],
                    carry_bits=1, rounded=True, booth=True, reuse_normalized_top=True,
                    output_extra_bits=oeb, defer_workspace=True, qubit_limit=qlim,
                    curve_keep_bits=4, defer_curve_phase=True, product_limit=plim,
                    borrow_integration=True
                )
                raw = Path(f"/tmp/test_var_{oeb}_{plim}_{qlim}.raw.kmx")
                compact = Path(f"/tmp/test_var_{oeb}_{plim}_{qlim}.compact.kmx")
                c.to_file(str(raw), format="kmx")
                subprocess.run([str(sys_path / 'compact_wires'), str(raw), str(compact), '22'], check=True, capture_output=True)
                c_comp = km.Circuit.from_file(str(compact))
                q = c_comp.num_qubits
                t = c_comp.operation_statistics().num_toffoli_like
                est_cost = q * (t + 250)
                print(f"oeb={oeb}, plim={plim}, qlim={qlim} -> Q={q}, T={t}, est_cost={est_cost:,}")
            except Exception as e:
                pass
