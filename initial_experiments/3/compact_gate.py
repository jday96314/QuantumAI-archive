"""Erase a gated inverse through its subfield norm, with less live workspace."""
import kickmix as km
from tower_inverse import tower,apply_map,linear_xor,product,del_inverse


class CompactGate:
    def _clear_inverse_scaffold(self,b,source,sc):
        info=tower(12)
        norm,invnorm,inner=sc[:6],sc[6:12],sc[12:]
        apply_map(b,source,info['ops'],True)
        a,c=source[:6],source[6:]
        del_inverse(b,norm,invnorm,inner)
        b.free(invnorm);b.free(inner)
        linear_xor(b,a,norm,info['square']);linear_xor(b,c,norm,info['musquare'])
        b.del_gf2_mul(a,c,target=norm,field=km.GF2Field(6));b.free(norm)
        apply_map(b,source,info['ops'])

    def _inverse_phase(self,b,source,mask,active=None):
        info=tower(12);field=km.GF2Field(6)
        apply_map(b,source,info['ops'],True)
        a,c=source[:6],source[6:]
        norm=b.alloc_qubits(6);product(b,a,c,norm)
        linear_xor(b,a,norm,info['square']);linear_xor(b,c,norm,info['musquare'])
        inv,sc=b.init_gf2_inverse_with_scaffold(norm,field=field)
        den=inv
        if active is not None:
            den=b.alloc_qubits(6);b.init_and(active,inv,den)
        b.cx(c,a);b.gf2_phase_by_product(a,den,mask=mask[:6],field=field);b.cx(c,a)
        b.gf2_phase_by_product(c,den,mask=mask[6:],field=field)
        if active is not None:
            b.del_and(active,inv,den);b.free(den)
        b.del_gf2_inverse_with_scaffold(norm,target=inv,scaffold=sc,field=field)
        b.free(inv);b.free(sc)
        linear_xor(b,a,norm,info['square']);linear_xor(b,c,norm,info['musquare'])
        b.del_gf2_mul(a,c,target=norm,field=field);b.free(norm)
        apply_map(b,source,info['ops'])

    def _gate(self,b,src,active):
        inv,sc=b.init_gf2_inverse_with_scaffold(src,field=self.gf)
        self._clear_inverse_scaffold(b,src,sc)
        g=b.alloc_qubits(12);b.init_and(active,inv,g)
        apply_map(b,inv,tower(12)['ops'],True)
        m=b.alloc_bits(12);b.hmr(inv,m);b.free(inv)
        self._inverse_phase(b,src,m)
        b.bit_store(m,False);b.free(m)
        return g

    def _ungate(self,b,src,active,g):
        apply_map(b,g,tower(12)['ops'],True)
        m=b.alloc_bits(12);b.hmr(g,m);b.free(g)
        self._inverse_phase(b,src,m,active)
        b.bit_store(m,False);b.free(m)


def test():
    import numpy as np
    from fast_dqi import FastBuilder
    class Gate(CompactGate):gf=km.GF2Field(12)
    gate=Gate();b=FastBuilder(km.CircuitBuilder())
    src=b.create_quantum_register(12,'src');active=b.create_quantum_register(1,'active')[0]
    out=b.create_quantum_register(12,'out')
    g=gate._gate(b,src,active);b.cx(g,out);gate._ungate(b,src,active,g)
    c=b.finish_circuit();N=8192
    values=np.arange(N)%4096;av=np.arange(N)//4096
    bits=lambda a,n:((a[None,:]>>np.arange(n)[:,None])&1).astype(bool)
    sim=km.Simulator(batch_size=N,seed=52);sim.use_same_registers_as(c)
    sim.write_across_shots(src,bits(values,12));sim.write_across_shots(km.array([active]),bits(av,1));sim.do(c)
    expected=np.array([gate.gf.pow(int(x),4094)*int(a) for x,a in zip(values,av)])
    assert np.array_equal(sim.read_across_shots(out),bits(expected,12))
    assert np.array_equal(sim.read_across_shots(src),bits(values,12))
    assert np.array_equal(sim.read_across_shots(km.array([active])),bits(av,1))
    public=set(list(src)+[active]+list(out))
    junk=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    assert not sim.read_across_shots(junk).any()
    assert not sim.read_phase_flipped_across_shots().any()
    print('PASS all inverse inputs and both activity values;',c.num_qubits,'qubits;',c.operation_statistics().num_toffoli_like,'Toffolis',flush=True)


if __name__=='__main__':test()
