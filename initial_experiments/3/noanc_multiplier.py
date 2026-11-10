"""Bilinear XOR multiplication using input/output conjugations and no scratch."""
import kickmix as km
from gf4_multiplier import terms
from crt_multiplier import indices


def apply_terms(b,lhs,rhs,target,triples):
    assert not (set(lhs)&set(rhs) or set(lhs)&set(target) or set(rhs)&set(target))
    for lm,rm,om in triples:
        li=list(indices(lm));ri=list(indices(rm));oi=list(indices(om))
        l,r,out=lhs[li[0]],rhs[ri[0]],target[oi[0]]
        ll=km.array([lhs[j] for j in li[1:]])
        rr=km.array([rhs[j] for j in ri[1:]])
        oo=km.array([target[j] for j in oi[1:]])
        b.cx(ll,l);b.cx(rr,r);b.cx(out,oo)
        b.ccx(l,r,out)
        b.cx(out,oo);b.cx(rr,r);b.cx(ll,l)


class NoAncMultiply:
    def _mul(self,b,x,y,target):
        apply_terms(b,x,y,target,terms(len(target)))


def test():
    import numpy as np
    for n in (6,12):
        f=km.GF2Field(n);rng=np.random.default_rng(n)
        b=km.CircuitBuilder();x=b.create_quantum_register(n,'x');y=b.create_quantum_register(n,'y');z=b.create_quantum_register(n,'z')
        apply_terms(b,x,y,z,terms(n));c=b.finish_circuit()
        N=4096;xx,yy,zz=rng.integers(0,1<<n,(3,N))
        for i in range(n):
            for j in range(n):xx[n*i+j],yy[n*i+j]=1<<i,1<<j
        bits=lambda a:((a[None,:]>>np.arange(n)[:,None])&1).astype(bool)
        sim=km.Simulator(batch_size=N,seed=93);sim.use_same_registers_as(c)
        for r,v in [(x,xx),(y,yy),(z,zz)]:sim.write_across_shots(r,bits(v))
        sim.do(c)
        expected=zz^np.array([f.mul(int(a),int(b)) for a,b in zip(xx,yy)])
        assert np.array_equal(sim.read_across_shots(z),bits(expected))
        assert np.array_equal(sim.read_across_shots(x),bits(xx))
        assert np.array_equal(sim.read_across_shots(y),bits(yy))
        assert not sim.read_phase_flipped_across_shots().any()
        assert c.num_qubits==3*n
        print('PASS no-scratch field product',n,c.operation_statistics().num_toffoli_like,flush=True)


if __name__=='__main__':test()
