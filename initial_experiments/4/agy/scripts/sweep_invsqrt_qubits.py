"""Sweep target_qubits in invsqrt_opt.py to find optimal score."""

import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers/invsqrt")
import kickmix as km
from invsqrt_opt import build_opt_invsqrt_circuit
from verify_submission_invsqrt import rotation_t_count, T_PER_TOFFOLI

def evaluate_invsqrt(circuit):
    stats = circuit.operation_statistics()
    toffoli_count = stats.num_toffoli_like
    t_base = 1.0
    n_rot = 35
    eps = 3.431142267989554e-13
    t_count = t_base + (n_rot * rotation_t_count(eps) if n_rot else 0)
    toffolis = toffoli_count + t_count / T_PER_TOFFOLI
    score = circuit.num_qubits * toffolis
    return circuit.num_qubits, toffolis, score

def main():
    base_score = 2076774.2337
    print(f"Baseline raw score: {base_score:,.2f}")
    for tq in [90, 85, 80, 75, 70, 68, 65, 62, 60, 58, 55, 52, 50]:
        try:
            c = build_opt_invsqrt_circuit(target_qubits=tq)
            q, t_eq, score = evaluate_invsqrt(c)
            ratio = base_score / score
            print(f"target_qubits={tq:2d} -> actual Q={q:2d}, T_eq={t_eq:.2f}, Score={score:,.2f}, Ratio={ratio:.4f}x")
        except Exception as e:
            print(f"target_qubits={tq:2d} -> Failed: {e}")

if __name__ == "__main__":
    main()
