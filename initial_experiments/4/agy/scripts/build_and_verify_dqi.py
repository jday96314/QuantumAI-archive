"""Build new dqi.kmb with n_iters=216 and verify with 1,536 shots."""

import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers/dqi")
import json
import time
import kickmix as km
from verify_submission_dqi import verify, ProblemSpecDQI
from dqi_opt import build_opt_dqi_circuit

def main():
    t0 = time.time()
    print("Building quantum circuit with default_n_iters=216...")
    c = build_opt_dqi_circuit()
    output_path = "/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/dqi.kmb"
    c.to_file(output_path, format="kmb")
    stats = c.operation_statistics()
    print(f"Wrote {output_path}: {c.num_qubits} qubits, {stats.num_toffoli_like:,} max Toffolis.")

    print("\nRunning verification on 1,536 shots with 4 workers...")
    spec = ProblemSpecDQI.benchmark_instance(num_shots=1536)
    res = verify(c, spec, max_workers=4)
    print("Verification passed!")
    for k, v in res.items():
        if k != "circuit":
            print(f"  {k}: {v}")

    base_dqi = 24439463583672.93
    ratio = base_dqi / res["score"]
    print(f"\nDQI Improvement Ratio: {ratio:.4f}x (Cost: {res['score']:,.2f})")

if __name__ == "__main__":
    main()
