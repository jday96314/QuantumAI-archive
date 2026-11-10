import kickmix as km
from fast_dqi import FastBuilder
from compact_dqi import CompactDecoder as Decoder
from dqi.dqi_reference_circuit import ProosAndZalkaDecoder as P
def build(iters,cls=Decoder):
    b=FastBuilder(km.CircuitBuilder())
    syn=[b.create_quantum_register(12,f'S{i}') for i in range(70)]
    mb=[b.create_classical_register(12,f'E{i}') for i in range(4095)]
    cls(70,km.GF2Field(12)).build_syndrome_decoding_mbuc(b,syndromes=syn,error_mbuc=mb,n_iters=iters)
    return b.finish_circuit().operation_statistics().num_toffoli_like
base=build(100)-build(99)
print('iter ~100 cost',base)
names=['_update_flag','_multiply_and_add_or_sub','_shift_b_array_right','_update_metadata','_swap_a_and_b','_update_halting_counter','_chien_search_and_forney_mbuc','_initial_shift_loop']
for nm in names:
    raw=P.__dict__[nm]
    setattr(P,nm,staticmethod(lambda *a,**k:None) if isinstance(raw,staticmethod) else (lambda *a,**k:None))
    if nm.startswith('_chien') or nm.startswith('_initial'):
        print(nm,'cost',build(1)-0,flush=True)
    else:
        print(nm,'per-iter',base-(build(100)-build(99)),flush=True)
    setattr(P,nm,raw)
print('build1',build(1))
