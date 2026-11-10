"""Clean, optimized linear inverse square root using Kickmix's native init_lookup and del_lookup."""

import sys
import time
from pathlib import Path
import numpy as np

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
    str(ROOT / "initial_experiments/1"),
]
import kickmix as km
from invsqrt import nr_reference_circuit as nr
import binary_centered_invsqrt as poly
from invsqrt.verify_submission_invsqrt import verify, ProblemSpecInvSqrt

def float_table(frac, order):
    h = 2.0**-frac
    table = np.zeros((2**(frac + 1), order + 1), dtype=float)
    for i in range(len(table)):
        a = 1 + (i % 2**frac + 0.5) * h
        scale = 1 if i >> frac else 1 / np.sqrt(2)
        cheb = np.polynomial.Chebyshev.interpolate(lambda x: scale / np.sqrt(a + x), order, domain=[-h/2, h/2])
        coef = cheb.convert(kind=np.polynomial.Polynomial).coef
        for j in range(order + 1):
            table[i, j] = coef[j] * (-1)**j
    return table

def table_for(frac, widths, bias):
    table = np.rint(float_table(frac, len(widths) - 1) * np.exp2(widths)).astype(np.int64)
    table[:, 0] += bias
    assert (table >= 0).all()
    return table

def build_fast_linear(frac=10, widths=(28, 18), bias=0):
    b = km.CircuitBuilder()
    x = b.create_quantum_register(22, 'x')
    
    # Normalization flags
    flags = b.alloc_qubits(10)
    nr._normalize(b, False, x, flags)
    
    s = b.alloc_qubits(4)
    for i in range(10):
        b.cx(flags[i], s[(9 - i) // 2])
    
    # Table data
    tab = table_for(frac, widths, bias)
    load_tab = tab[:, [0, 1]]
    
    target_c0 = b.alloc_qubits(widths[0])
    target_c1 = b.alloc_qubits(widths[1])
    target_qubits = km.array([*target_c0, *target_c1])
    
    # Address: top frac bits of normalized x plus parity bit s[0]
    address = km.array([*x[21 - frac:21], flags[0]])
    
    # Pack table for kickmix lookup
    packed = [int(row[0]) | (int(row[1]) << widths[0]) for row in load_tab]
    
    # 1. Native fast init_lookup
    b.init_lookup(table=packed, address=address, target=target_qubits)
    
    # 2. Subtract c1 * Delta x from c0
    # In binary-centered normalization:
    # Delta x is centered: z = (x & ((1<<L)-1)) ^ (1<<(L-1))
    L = 21 - frac
    z = b.alloc_qubits(L)
    b.cx(x[:L], z)
    b.x(z[-1]) # centered flip
    
    # Booth multiplication of c1 by z into target_c0
    # Or linear addition
    # Here Delta x has only 21 - 10 = 11 bits!
    # c1 has 18 bits.
    from booth_adder import accumulate
    accumulate(
        b, target_c1, z, target_c0,
        widths[0] - widths[1] - 21,
        subtract=True, rounded=True, booth=True
    )
    
    # Controlled shift by s
    from binary_centered_invsqrt import cshift
    y = target_c0
    for j in range(4):
        cshift(b, y, 1 << j, s[j], len(y))
    
    # Apply phase via Z and Z_POW
    nr._phase_by_bits(b, y[:widths[0] + 10])
    
    # Reverse shift
    for j in reversed(range(4)):
        cshift(b, y, 1 << j, s[j], len(y), inv=True)
    
    # Reverse Booth multiplication
    accumulate(
        b, target_c1, z, target_c0,
        widths[0] - widths[1] - 21,
        subtract=False, rounded=True, booth=True
    )
    b.x(z[-1])
    b.cx(x[:L], z)
    b.free(z)
    
    # 3. Native fast del_lookup
    b.del_lookup(table=packed, address=address, target=target_qubits)
    b.free(target_qubits)
    
    for i in reversed(range(10)):
        b.cx(flags[i], s[(9 - i) // 2])
    b.free(s)
    
    nr._normalize(b, True, x, flags)
    b.free(flags)
    
    c = b.finish_circuit()
    return c

print("Building fast linear circuit...")
start = time.time()
c = build_fast_linear(10, (28, 18), 0)
print(f"Built in {time.time()-start:.2f}s!")
print(f"Total qubits: {c.num_qubits}")
print(f"Toffoli-like: {c.operation_statistics().num_toffoli_like}")

out_path = Path("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/runs/fast_linear_invsqrt.kmb")
c.to_file(str(out_path), format="kmb")
print(f"Saved to {out_path}")
