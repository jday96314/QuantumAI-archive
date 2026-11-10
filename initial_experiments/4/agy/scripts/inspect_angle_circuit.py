"""Inspect operations and structure in baseline angle_k_0.kmb."""

import kickmix as km
import sys
sys.path.insert(0, "/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers")
from dfthcsos.scoreinfra.circuits import register_layout

def main():
    c = km.Circuit.from_file("/media/shared/ExpansionDrive2/QuantumAI/from_host/sample_submission/dfthcsos/angle_k_0.kmb")
    print("Num qubits:", c.num_qubits)
    layout = register_layout([c])
    print("Selection:", len(layout.selection), "qubits")
    print("Words:", len(layout.words), "registers of", layout.brot, "qubits")
    stats = c.operation_statistics()
    print("Num toffoli like:", stats.num_toffoli_like)
    print("Num ops:", len(c.operations))
    
    # Check first 20 operations
    for i, op in enumerate(c.operations[:20]):
        print(f"op[{i}]:", op)

if __name__ == "__main__":
    main()
