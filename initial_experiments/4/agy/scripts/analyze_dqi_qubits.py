"""Profile exact qubit allocations and peak concurrency in dqi_opt.py."""

import kickmix as km
from dqi_opt import build_opt_dqi_circuit

def main():
    circuit = build_opt_dqi_circuit(n_iters=210)
    print("Num qubits:", circuit.num_qubits)
    stats = circuit.operation_statistics()
    print("Num toffolis:", stats.num_toffoli_like)

if __name__ == "__main__":
    main()
