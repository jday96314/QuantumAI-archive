"""Verify the full quantum circuit with DQI early halting and shift_b_array_right optimization."""

from __future__ import annotations

import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers/dqi")
import kickmix as km
import numpy as np
from verify_submission_dqi import verify, ProblemSpecDQI
from dqi_opt import build_opt_dqi_circuit

def main():
    print("Building quantum circuit with default_n_iters=216 and optimized shift...")
    c = build_opt_dqi_circuit()
    print(f"Built circuit: {c.num_qubits} qubits, {c.operation_statistics().num_toffoli_like:,} max Toffolis.")

    print("\nRunning verification on 256 shots with 4 workers...")
    spec = ProblemSpecDQI.benchmark_instance(num_shots=256)
    res = verify(c, spec, max_workers=4)
    print("Verification result:")
    for k, v in res.items():
        print(f"  {k}: {v}")

if __name__ == "__main__":
    main()
