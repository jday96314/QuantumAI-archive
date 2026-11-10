"""Test replacing unary_iteration with b.init_lookup in binary_centered_invsqrt."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
    str(ROOT / "initial_experiments/1"),
]
import kickmix as km
import numpy as np

# Copy binary_centered_invsqrt logic but replace lookup with init_lookup
import binary_centered_invsqrt as poly

# Let's inspect poly.build source or create a patched build function
import inspect
src = inspect.getsource(poly.build)

# Replace unary_iteration block with init_lookup
old_block = """        with b.unary_iteration(address) as cursor:
            for a in range(len(tab)):
                cursor.move_to(a)
                for j,r in enumerate(load_regs):b.gf2_iadd(int(load_tab[a,j]),target=r,control=cursor.match_qubit)"""

new_block = """        load_sizes = [len(r) for r in load_regs]
        packed = [sum(int(v) << sum(load_sizes[:k]) for k,v in enumerate(row)) for row in load_tab]
        target = km.array([q for reg in load_regs for q in reg])
        b.init_lookup(table=packed, address=address, target=target)"""

assert old_block in src, "old_block not found in poly.build"
new_src = src.replace(old_block, new_block)

# Execute new build function in poly's module dict
mod_dict = dict(poly.__dict__)
exec(new_src, mod_dict)
patched_build = mod_dict['build']

print("Building linear invsqrt with init_lookup...")
c = patched_build(10, (28, 18), bias=0, rounded=True, booth=True)
print(f"Success! num_qubits={c.num_qubits}, toffoli_like={c.operation_statistics().num_toffoli_like}")

# Save and verify on small n_p
from invsqrt.verify_submission_invsqrt import verify, ProblemSpecInvSqrt
spec = ProblemSpecInvSqrt.for_n_p(6)
res = verify(c, spec, max_workers=4)
print("Verification on n_p=6:", res["status"], "score:", res.get("score"), "worst_ratio:", res.get("worst_error_over_tolerance"))
