"""Build and profile new DQI circuit with breakthrough logic."""

import kickmix as km
from dqi_opt import build_opt_dqi_circuit

def main():
    print("Building DQI circuit with n_iters=210...")
    c = build_opt_dqi_circuit(n_iters=210)
    print(f"Num qubits: {c.num_qubits}")
    stats = c.operation_statistics()
    print(f"Num Toffolis: {stats.num_toffoli_like:,}")
    
    # Calculate score assuming 0 errors (rep_factor = 1.19126)
    rep_factor = 1.1912648483863488
    # In full DQI, toffolis per shot is the average
    # Baseline: base_toff = 5725108.38, base_score = 24439463583672.93
    base_score = 24439463583672.93
    # With n_iters=210, let's see the ratio
    max_toff = stats.num_toffoli_like
    print(f"Max Toffolis: {max_toff:,}")

if __name__ == "__main__":
    main()
