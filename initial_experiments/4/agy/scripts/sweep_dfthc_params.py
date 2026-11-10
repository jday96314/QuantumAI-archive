"""Search parameter space of dfthcsos for minimum theoretical cost."""

import math
import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers")
from dfthcsos.verify import heuristic_costs

N = 76

def evaluate_cost(R, B, C, brot=12, J=8):
    A = N + R * B
    # Selection register bits
    sel = math.ceil(math.log2(A))
    # Angle circuits count
    K = math.ceil((N - 1) / J)
    # Check qubit cap
    ang_qubits = sel + J * brot + (sel - 1)
    if ang_qubits > 128:
        return None
    
    # Angle Toffolis (ordinary QROM per circuit, 2 passes)
    ang_toff = 2 * K * (A - 2)
    
    try:
        h = heuristic_costs(N, R, B, C, brot)
    except Exception:
        return None
        
    T_total = ang_toff + h["toffoli"]
    Q_total = h["persistent_qubits"] + (J * brot) + max(h["temporary_qubits"], sel - 1)
    
    return T_total, Q_total, T_total * Q_total, h

def main():
    base_T = 93237
    base_Q = 730
    base_cost = base_T * base_Q * 49.30136925057439 # 3,355,599,588.32
    print(f"Baseline: T={base_T}, Q={base_Q}, Cost={base_cost:,.2f}")
    
    # Baseline shape: R=65, B=60, C=60
    # What if we change C? (Remember: C can be compressed via SVD without changing U or B!)
    print("\n--- Sweeping C with fixed R=65, B=60, brot=12, J=8 ---")
    for C in [60, 50, 40, 30, 25, 20, 16, 12, 10, 8, 4, 2, 1]:
        res = evaluate_cost(65, 60, C, brot=12, J=8)
        if res:
            T, Q, TQ, h = res
            cost = 49.3 * TQ
            ratio = base_cost / cost
            print(f"C={C:2d}: T={T:6d} (ang={2*10*(3974):5d}, heur={h['toffoli']:5d}), Q={Q:3d} (temp={h['temporary_qubits']:3d}), Cost={cost:,.0f}, Ratio={ratio:.3f}x")

    print("\n--- Sweeping B with fixed R=65, C=20, brot=12, J=8 ---")
    for B in [60, 50, 40, 32, 24, 20, 16, 12, 8]:
        res = evaluate_cost(65, B, 20, brot=12, J=8)
        if res:
            T, Q, TQ, h = res
            cost = 49.3 * TQ
            ratio = base_cost / cost
            print(f"B={B:2d}: T={T:6d}, Q={Q:3d}, Cost={cost:,.0f}, Ratio={ratio:.3f}x")

    print("\n--- Sweeping R with B=16, C=10, brot=12, J=8 ---")
    for R in [65, 50, 40, 32, 25, 20, 16, 10, 5, 2, 1]:
        res = evaluate_cost(R, 16, 10, brot=12, J=8)
        if res:
            T, Q, TQ, h = res
            cost = 49.3 * TQ
            ratio = base_cost / cost
            print(f"R={R:2d}: T={T:6d}, Q={Q:3d}, Cost={cost:,.0f}, Ratio={ratio:.3f}x")

if __name__ == "__main__":
    main()
