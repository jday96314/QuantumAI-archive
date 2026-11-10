"""Compute one controlled head correction and reuse it over an extraction row."""
from functools import lru_cache
import kickmix as km
from extension_inverse import Adapter,setup
from schur_control import caxpy_rt
from fast_dqi import FastBuilder

@lru_cache(None)
def build(constants,erase=False):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());flags=b.create_quantum_register(len(constants),'flags');src=b.create_quantum_register(12,'src');dst=b.create_quantum_register(12,'dst');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));rt=[]
 for flag,c in zip(flags,constants):rt+=caxpy_rt(m,flag.id,[q.id for q in src],[q.id for q in dst],c)
 (E.free if erase else E.mac)(rt);ab.release_recycle('shared_head_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))
