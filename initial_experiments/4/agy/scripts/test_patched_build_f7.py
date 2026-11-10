"""Test replacing unary_iteration with b.init_lookup in binary_centered_invsqrt on frac=7."""

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
import inspect
import binary_centered_invsqrt as poly

# Load candidate table
table = np.load(sys_path / 'candidates/margin_curve8_t5.9e-08.npy')
poly.table_for = lambda *args: table

src = inspect.getsource(poly.build)

old_block = """        with b.unary_iteration(address) as cursor:
            for a in range(len(tab)):
                cursor.move_to(a)
                for j,r in enumerate(load_regs):b.gf2_iadd(int(load_tab[a,j]),target=r,control=cursor.match_qubit)"""

new_block = """        load_sizes = [len(r) for r in load_regs]
        packed = [sum(int(v) << sum(load_sizes[:k]) for k,v in enumerate(row)) for row in load_tab]
        target = km.array([q for reg in load_regs for q in reg])
        b.init_lookup(table=packed, address=address, target=target)"""

assert old_block in src
new_src = src.replace(old_block, new_block)

mod_dict = dict(poly.__dict__)
exec(new_src, mod_dict)
patched_build = mod_dict['build']

print("Building frac=7 with init_lookup...")
c = patched_build(
    7, (27, 18, 8), lookup_frac=8, accum_bits=[18, 9, None],
    carry_bits=1, rounded=True, booth=True, reuse_normalized_top=True,
    output_extra_bits=10, defer_workspace=True, qubit_limit=94,
    curve_keep_bits=4, defer_curve_phase=True, product_limit=17,
    borrow_integration=True
)

print(f"Patched build: qubits={c.num_qubits}, toffoli={c.operation_statistics().num_toffoli_like}")

# Now compact wires
raw = Path("/tmp/patched_f7.raw.kmx")
compact = Path("/tmp/patched_f7.compact.kmx")
c.to_file(str(raw), format="kmx")
import subprocess
subprocess.run([str(sys_path / 'compact_wires'), str(raw), str(compact), '22'], check=True)
c_comp = km.Circuit.from_file(str(compact))
print(f"Compacted build: qubits={c_comp.num_qubits}, toffoli={c_comp.operation_statistics().num_toffoli_like}")
