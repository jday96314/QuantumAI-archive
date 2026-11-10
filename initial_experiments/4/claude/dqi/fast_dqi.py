"""Production decoder builder using GF(4) multiplication and subfield inversion."""
import kickmix as km
import crt_multiplier as crt
import tower_inverse as tower
from gf4_multiplier import terms as gf4_terms

# The original construction functions remain available for independent algebra checks.
# These providers select exact binary linear forms; all emitted instructions are kickmix gates.
legacy_tower_terms=tower.crt_terms

def selected_terms(n):
    return gf4_terms(n) if n in [6,12] else legacy_tower_terms(n)


def install():
    crt.terms12=lambda modulus=None:gf4_terms(12)
    tower.crt_terms=selected_terms


class FastBuilder(tower.TowerBuilder):
    def __init__(self,builder):
        install()
        super().__init__(builder)


def test():
    install()
    tower.test()
    crt.test()
    # Division integration and all input branches of the final builder.
    import numpy as np
    rng=np.random.default_rng(31)
    b=FastBuilder(km.CircuitBuilder());a=b.create_quantum_register(12,'a');v=b.create_quantum_register(12,'v')
    q,sc=b.init_gf2_div_with_scaffold(a,v)
    out=b.create_quantum_register(12,'out');b.cx(q,out)
    b.del_gf2_div_with_scaffold(a,v,target=q,scaffold=sc);b.free(q);b.free(sc)
    c=b.finish_circuit();N=4096;s=km.Simulator(batch_size=N,seed=41);s.use_same_registers_as(c)
    aa=rng.integers(0,4096,N);vv=np.arange(N)
    bits=lambda x:((x[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    s.write_across_shots(a,bits(aa));s.write_across_shots(v,bits(vv));s.do(c)
    f=km.GF2Field(12);expected=np.array([f.mul(int(x),f.pow(int(y),4094)) for x,y in zip(aa,vv)])
    assert np.array_equal(s.read_across_shots(out),bits(expected))
    assert np.array_equal(s.read_across_shots(a),bits(aa))
    assert np.array_equal(s.read_across_shots(v),bits(vv))
    assert not s.read_phase_flipped_across_shots().any()
    public={q.id for r in [a,v,out] for q in r}
    assert not s.read_across_shots(km.array([km.q(i) for i in range(c.num_qubits) if i not in public])).any()
    print('Division PASS every divisor, Toffolis',c.operation_statistics().num_toffoli_like,flush=True)

if __name__=='__main__':test()
