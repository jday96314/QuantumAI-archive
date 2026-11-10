"""Tabulate T_total / Q_total of the dfthcsos cost model across factorization shapes.

Angle-lookup cost is modelled as a plain unary-iteration QROM per batch
(A - 2 Toffolis per circuit, two passes), with J the largest batch that fits
in 128 qubits (selection + J*brot + ~selection-1 unary ancillas).
"""

from __future__ import annotations

import math
import sys

sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers")
from dfthcsos.verify import heuristic_costs  # noqa: E402

N = 76


def angle_cost(A: int, brot: int, J: int) -> tuple[int, int, int]:
    sel = math.ceil(math.log2(A))
    K = math.ceil((N - 1) / J)
    tof = 2 * K * (A - 2)
    qubits = sel + J * brot + sel
    return tof, J * brot, qubits - sel - J * brot


def best(R: int, B: int, C: int, brot: int, lam: float = 49.3) -> tuple:
    A = N + R * B
    h = heuristic_costs(N, R, B, C, brot)
    out = None
    sel = math.ceil(math.log2(A))
    for J in range(1, 76):
        if sel + J * brot + sel > 128:
            break
        at, ap, atemp = angle_cost(A, brot, J)
        T = at + h["toffoli"]
        Q = h["persistent_qubits"] + ap + max(h["temporary_qubits"], atemp)
        s = lam * T * Q
        if out is None or s < out[0]:
            out = (s, J, T, Q, at, h["toffoli"], h["persistent_qubits"], h["temporary_qubits"])
    return out


def main() -> None:
    base = 3355599588.3155375
    shapes = [(65, 60, 60), (65, 60, 20), (65, 40, 40), (1, 300, 300), (1, 400, 400), (1, 500, 300),
              (4, 150, 100), (10, 100, 50), (20, 60, 30), (2, 250, 250)]
    for brot in (12, 11, 10):
        for R, B, C in shapes:
            s, J, T, Q, at, ht, hp, ht2 = best(R, B, C, brot)
            print(f"brot={brot} R={R:3d} B={B:3d} C={C:3d}: J={J:2d} T={T:7d} (ang {at:6d}, heur {ht:6d}) "
                  f"Q={Q:4d} (pers {hp}, temp {ht2})  ratio@lam49.3={base / s:6.2f}")


if __name__ == "__main__":
    main()
