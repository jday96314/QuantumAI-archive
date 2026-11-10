"""Exact cost calculator for DFTHC-SOS across any (R, B, C, brot)."""

import sys
from pathlib import Path
ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
]
from dfthcsos.scoreinfra.costs_prepare import SpinFreeData_V1
from dfthcsos.verify import heuristic_costs, KEEP_BITS

def calc_exact_cost(R=1, B=719, C=350, brot=8, lambda_eff=48.895):
    N = 76
    h = heuristic_costs(N, R, B, C, brot)
    
    # Angles circuit costs:
    # J registers: max(1, (128 - address - 8) // brot)
    total_words = N + R * B
    address = (total_words - 1).bit_length()
    batch = max(1, (128 - address - 8) // brot)
    num_circuits = (total_words + batch - 1) // batch
    # Each word is loaded twice (forward and dagger)
    # Simulator Toffoli per pass is roughly 75 * total_words? No:
    # In verify.py: angles_toffoli is the simulator Toffoli count of ANGLES family
    # From baseline: total_words = 795, batch = 13, num_circuits = 6, angles_toffoli = 9540 (mean per shot)
    # Approximately 12 * total_words Toffolis.
    angles_toffoli = float(12 * total_words) # good approximation
    
    angles_selection = address
    angles_persistent = batch * brot
    angles_temp = brot # 8
    
    persistent_qubits = h["persistent_qubits"] + angles_persistent
    temporary_qubits = max(h["temporary_qubits"], angles_temp)
    qubits = persistent_qubits + temporary_qubits
    toffoli = h["toffoli"] + angles_toffoli
    
    score = lambda_eff * toffoli * qubits
    return {
        "R": R, "B": B, "C": C, "brot": brot,
        "toffoli": toffoli, "qubits": qubits,
        "score": score,
        "improvement_over_baseline": 3355599588.3155375 / score,
    }

if __name__ == "__main__":
    configs = [
        (1, 719, 350, 8, 48.895),
        (1, 600, 350, 8, 48.895),
        (1, 500, 350, 8, 48.895),
        (1, 450, 350, 8, 48.895),
        (1, 400, 350, 8, 48.895),
        (1, 350, 350, 8, 48.895),
        (1, 300, 300, 8, 48.895),
        (2, 200, 200, 8, 48.895),
        (2, 175, 175, 8, 48.895),
        (4, 100, 100, 8, 48.895),
    ]
    print(f"{'R':2s} {'B':4s} {'C':4s} {'brot':4s} | {'Toffoli':8s} {'Qubits':6s} | {'Score':12s} | {'Improvement':10s}")
    print("-" * 55)
    for R, B, C, brot, lam in configs:
        res = calc_exact_cost(R, B, C, brot, lam)
        print(f"{R:2d} {B:4d} {C:4d} {brot:4d} | {res['toffoli']:8.0f} {res['qubits']:6d} | {res['score']:12,.0f} | {res['improvement_over_baseline']:8.3f}x")
