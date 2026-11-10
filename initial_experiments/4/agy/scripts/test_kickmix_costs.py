"""Check what operations count as toffoli_like in kickmix."""

import kickmix as km

def main():
    b = km.CircuitBuilder()
    q = b.alloc_qubits(5)
    b.init_and(q[0], q[1], q[2])
    c1 = b.finish_circuit()
    print("init_and num_toffoli_like:", c1.operation_statistics().num_toffoli_like)

    b = km.CircuitBuilder()
    q = b.alloc_qubits(5)
    b.del_and(q[0], q[1], q[2])
    c2 = b.finish_circuit()
    print("del_and num_toffoli_like:", c2.operation_statistics().num_toffoli_like)

    b = km.CircuitBuilder()
    q = b.alloc_qubits(5)
    b.ccx(q[0], q[1], q[2])
    c3 = b.finish_circuit()
    print("ccx num_toffoli_like:", c3.operation_statistics().num_toffoli_like)

    b = km.CircuitBuilder()
    q = b.alloc_qubits(5)
    b.cx(q[0], q[1])
    c4 = b.finish_circuit()
    print("cx num_toffoli_like:", c4.operation_statistics().num_toffoli_like)

if __name__ == "__main__":
    main()
