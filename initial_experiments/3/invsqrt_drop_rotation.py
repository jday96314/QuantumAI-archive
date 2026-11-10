"""Build the fitted inverse-square-root circuit without its least phase bit."""
import argparse
from fractions import Fraction
from pathlib import Path

import numpy as np

import polynomial_invsqrt as poly
from invsqrt import nr_reference_circuit as nr


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    table = np.load('initial_experiments/3/parameters/invsqrt_table.npy')
    original_table = poly.table_for
    original_phase = nr._phase_by_bits

    def fitted(f, w, bias=0):
        return table if f == 7 and tuple(w) == (28, 19, 8) else original_table(f, w, bias)

    def phase_without_lowest(builder, y):
        f = len(y)
        for j in range(1, f):
            builder.z_pow(y[j], Fraction(2) ** (j - f + 1))

    poly.table_for = fitted
    nr._phase_by_bits = phase_without_lowest
    try:
        circuit = poly.build(7, (28, 19, 8), 0)
    finally:
        poly.table_for = original_table
        nr._phase_by_bits = original_phase
    a.output.parent.mkdir(parents=True, exist_ok=True)
    circuit.to_file(a.output, format='kmb')
    print(circuit.num_qubits, circuit.operation_statistics().num_toffoli_like)


if __name__ == '__main__':
    main()
