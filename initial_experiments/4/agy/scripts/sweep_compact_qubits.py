"""Sweep qubit_limit on build_compact_invsqrt to see if fewer qubits are possible."""

from pathlib import Path
import subprocess
import numpy as np
import kickmix as km

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys_path = ROOT / "initial_experiments/3"
import sys
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(sys_path),
    str(ROOT / "initial_experiments/1"),
]
import binary_centered_invsqrt as polynomial

table = np.load(sys_path / 'candidates/margin_curve8_t5.9e-08.npy')
polynomial.table_for = lambda *args: table

for qlim in [94, 92, 90, 88, 86, 84, 82, 80, 78, 76, 75]:
    try:
        c = polynomial.build(
            7, (27, 18, 8), lookup_frac=8, accum_bits=[18, 9, None],
            carry_bits=1, rounded=True, booth=True, reuse_normalized_top=True,
            output_extra_bits=10, defer_workspace=True, qubit_limit=qlim,
            curve_keep_bits=4, defer_curve_phase=True, product_limit=17,
            borrow_integration=True
        )
        # Compact wires
        raw = Path(f"/tmp/test_qlim_{qlim}.raw.kmx")
        compact = Path(f"/tmp/test_qlim_{qlim}.compact.kmx")
        c.to_file(str(raw), format="kmx")
        res = subprocess.run([str(sys_path / 'compact_wires'), str(raw), str(compact), '22'], capture_output=True, text=True)
        c_comp = km.Circuit.from_file(str(compact))
        q = c_comp.num_qubits
        t = c_comp.operation_statistics().num_toffoli_like
        score = q * (t + 250)
        print(f"qubit_limit={qlim:2d} -> alloc={c.num_qubits}, compacted={q}, Toffoli={t}, Est Score={score:,} | {res.stdout.strip()}")
    except Exception as e:
        print(f"qubit_limit={qlim:2d} -> FAILED: {e}")
