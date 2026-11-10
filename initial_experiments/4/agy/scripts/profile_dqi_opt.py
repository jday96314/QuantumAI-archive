"""Profile Toffoli and qubit breakdown of dqi.kmb by mark."""

import kickmix as km


def main() -> None:
    c = km.Circuit.from_file("/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/4/agy/dqi.kmb")
    stats = c.operation_statistics()
    print(f"Total qubits: {c.num_qubits}")
    print(f"Total Toffoli-like: {stats.num_toffoli_like:,}")
    for name, st in sorted(stats.items(), key=lambda kv: -kv[1].num_toffoli_like):
        print(f"  {str(name):35s}: {st.num_toffoli_like:12,d} Toffolis")


if __name__ == "__main__":
    main()
