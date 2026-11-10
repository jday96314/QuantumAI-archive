import sys,kickmix as km
from fast_dqi import FastBuilder
from compact_dqi import CompactDecoder as Decoder
def build(iters):
    b=FastBuilder(km.CircuitBuilder())
    syn=[b.create_quantum_register(12,f'S{i}') for i in range(70)]
    mb=[b.create_classical_register(12,f'E{i}') for i in range(4095)]
    Decoder(70,km.GF2Field(12)).build_syndrome_decoding_mbuc(b,syndromes=syn,error_mbuc=mb,n_iters=iters)
    c=b.finish_circuit();return c
for it in [0,1,2]:
    c=build(it);print(it,c.operation_statistics().num_toffoli_like,c.num_qubits,flush=True)
