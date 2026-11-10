"""Search extreme parameter space for DFTHC-SOS."""

import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers")
from dfthcsos.verify import heuristic_costs
import math

N = 76
base_cost = 3355599588.3155375

def eval_config(R, B, C, brot=12, J=8):
    A = N + R * B
    sel = math.ceil(math.log2(A))
    K = math.ceil((N - 1) / J)
    ang_qubits = sel + J * brot + (sel - 1)
    if ang_qubits > 128:
        return None
    ang_toff = 2 * K * (A - 2)
    try:
        h = heuristic_costs(N, R, B, C, brot)
    except Exception:
        return None
    T_total = ang_toff + h["toffoli"]
    Q_total = h["persistent_qubits"] + (J * brot) + max(h["temporary_qubits"], sel - 1)
    cost = 49.30136925057439 * T_total * Q_total
    ratio = base_cost / cost
    return T_total, Q_total, cost, ratio

def main():
    print("--- Searching extreme configurations ---")
    best = []
    for R in [65, 50, 40, 30, 20, 15, 10, 8, 6, 5, 4, 3, 2, 1]:
        for B in [60, 40, 32, 24, 20, 16, 12, 8, 6, 4, 2, 1]:
            for C in [60, 40, 30, 20, 16, 12, 10, 8, 6, 4, 2, 1]:
                for J in [8, 9, 10, 11, 12, 13, 14, 15]:
                    res = eval_config(R, B, C, brot=12, J=J)
                    if res:
                        T, Q, cost, ratio = res
                        best.append((ratio, R, B, C, J, T, Q, cost))
    best.sort(key=lambda x: x[0], reverse=True)
    print("Top 20 configurations:")
    for ratio, R, B, C, J, T, Q, cost in best[:20]:
        print(f"R={R:2d}, B={B:2d}, C={C:2d}, J={J:2d} -> T={T:5d}, Q={Q:3d}, Cost={cost:,.0f}, Ratio={ratio:.2f}x")

if __name__ == "__main__":
    main()
