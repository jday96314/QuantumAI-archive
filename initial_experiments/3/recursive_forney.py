"""Contract inverse parity down the complete GF(4096)/GF(64)/GF(8) tower."""
from functools import lru_cache
import kickmix as km
from tower_inverse import tower,apply_map,map_ops,product,linear_xor
from forney_dqi import product_terms
from transposed_forney import TransposedForney,transpose_product,transpose_phase


@lru_cache(None)
def trans_ops(n):
    cols=tower(n)['columns']
    trans=tuple(sum(((cols[i]>>j)&1)<<i for i in range(n)) for j in range(n))
    return map_ops(trans)


def inverse_dot_data(b,source,u):
    n=len(source)
    if n==3:
        inv,sc=b.init_gf2_inverse_with_scaffold(source,field=km.GF2Field(3))
        def cleanup():
            b.del_gf2_inverse_with_scaffold(source,target=inv,scaffold=sc,field=km.GF2Field(3))
            b.free(inv);b.free(sc)
        return u,inv,cleanup
    info=tower(n);m=n//2;field=km.GF2Field(m)
    apply_map(b,source,info['ops'],True);apply_map(b,u,trans_ops(n))
    a,c=source[:m],source[m:]
    norm=b.alloc_qubits(m);v=b.alloc_qubits(m)
    product(b,a,c,norm)
    linear_xor(b,a,norm,info['square']);linear_xor(b,c,norm,info['musquare'])
    b.cx(c,a);transpose_product(b,a,u[:m],v);b.cx(c,a)
    transpose_product(b,c,u[m:],v)
    leaf,inv,inner_cleanup=inverse_dot_data(b,norm,v)
    def cleanup():
        inner_cleanup()
        vm=b.alloc_bits(m);b.hmr(v,vm)
        b.cx(c,a);transpose_phase(b,a,u[:m],vm);b.cx(c,a)
        transpose_phase(b,c,u[m:],vm)
        b.bit_store(vm,False);b.free(vm);b.free(v)
        linear_xor(b,a,norm,info['square']);linear_xor(b,c,norm,info['musquare'])
        b.del_gf2_mul(a,c,target=norm,field=field);b.free(norm)
        apply_map(b,u,trans_ops(n),True);apply_map(b,source,info['ops'])
    return leaf,inv,cleanup


class RecursiveForney(TransposedForney):
    def _forney_kernel(self,b,s,d,o,mask):
        tterms=product_terms(self.gf)
        u=b.alloc_qubits(12)
        def mask_map():
            for k in range(12):
                b.push_condition(mask[k])
                b.cx(km.array([o[i] for i,j in tterms[k]]),km.array([u[j] for i,j in tterms[k]]))
                b.pop_condition()
        mask_map()
        leaf,inv,cleanup=inverse_dot_data(b,d,u)
        parity=b.alloc_qubits(1)[0];ladder=b.alloc_qubits(10);meas=b.alloc_bits(1)[0]
        for j in range(3):b.ccx(leaf[j],inv[j],parity)
        b.x(s);b.init_and(s[0],s[1],ladder[0])
        for j in range(1,10):b.init_and(ladder[j-1],s[j+1],ladder[j])
        b.ccz(ladder[9],s[11],parity)
        for j in reversed(range(1,10)):b.del_and(ladder[j-1],s[j+1],ladder[j])
        b.del_and(s[0],s[1],ladder[0]);b.x(s)
        b.hmr(parity,meas);b.push_condition(meas);b.cz(leaf,inv);b.pop_condition()
        cleanup();mask_map()
        b.bit_store(meas,False);b.free(meas)
        b.free(u);b.free(ladder);b.free(parity)


def test():
    import numpy as np
    from fast_dqi import FastBuilder
    b=FastBuilder(km.CircuitBuilder())
    s=b.create_quantum_register(12,'s');d=b.create_quantum_register(12,'d')
    o=b.create_quantum_register(12,'o');mask=b.create_classical_register(12,'m')
    class Kernel(RecursiveForney):gf=km.GF2Field(12)
    Kernel()._forney_kernel(b,s,d,o,mask)
    c=b.finish_circuit();N=4096;rng=np.random.default_rng(284)
    ss=rng.integers(0,4096,N);ss[::2]=0
    dd=np.arange(N);oo=rng.integers(0,4096,N);mm=rng.integers(0,4096,N)
    bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    sim=km.Simulator(batch_size=N,seed=38);sim.use_same_registers_as(c)
    for reg,val in [(s,ss),(d,dd),(o,oo),(mask,mm)]:sim.write_across_shots(reg,bits(val))
    sim.do(c)
    gf=Kernel.gf
    expected=np.array([s0==0 and (gf.mul(int(o0),gf.pow(int(d0),4094))&int(m0)).bit_count()%2==1 for s0,d0,o0,m0 in zip(ss,dd,oo,mm)])
    assert np.array_equal(sim.read_phase_flipped_across_shots(),expected)
    for reg,val in [(s,ss),(d,dd),(o,oo)]:assert np.array_equal(sim.read_across_shots(reg),bits(val))
    public=set(list(s)+list(d)+list(o))
    junk=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    assert not sim.read_across_shots(junk).any()
    print('PASS recursive inverse parity all 4096 divisors;',c.num_qubits,'qubits;',c.operation_statistics().num_toffoli_like,'Toffolis',flush=True)


if __name__=='__main__':test()
