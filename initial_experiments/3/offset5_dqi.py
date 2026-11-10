"""Repaired EEA using no-scratch products and compact Forney arithmetic."""
import kickmix as km
from offset_dqi import OffsetDecoder
from transposed_forney import TransposedForney
from compact_gate import CompactGate
from noanc_multiplier import NoAncMultiply


class Decoder(NoAncMultiply,CompactGate,TransposedForney,OffsetDecoder):
    pass


def build_circuit(m=4095,n=70,chien=True):
    from fast_dqi import FastBuilder
    b=FastBuilder(km.CircuitBuilder())
    syn=[b.create_quantum_register(12,f'S{i}') for i in range(n)]
    mb=[b.create_classical_register(12,f'E{i}') for i in range(m)]
    Decoder(n,km.GF2Field(12)).build(b,syn,mb,chien)
    return b.finish_circuit()


if __name__=='__main__':
    import sys,time
    start=time.time();c=build_circuit()
    print('built',c.num_qubits,c.operation_statistics().num_toffoli_like,time.time()-start,flush=True)
    c.to_file(sys.argv[1],format='kmb')
