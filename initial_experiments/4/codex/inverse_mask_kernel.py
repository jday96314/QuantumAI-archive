"""Root-gated parity of a field inverse against a runtime classical mask.

A candidate primitive for a Wronskian-based DQI scan. This does not by itself
replace the decoder: the reversible in-place Wronskian transform is unfinished.
"""
import json
import kickmix as km
from tower_inverse import tower,apply_map,product,linear_xor
from recursive_forney import inverse_dot_data,trans_ops
from transposed_forney import multiplication_terms
from crt_multiplier import indices


def classical_transpose(b,x,u,target):
    flag=b.alloc_bits(1)[0]
    for lm,rm,om in multiplication_terms(len(x)):
        xx=list(indices(lm));lead=x[xx[0]];other=km.array([x[j] for j in xx[1:]])
        cm=km.array([u[j] for j in indices(om)]);dst=km.array([target[j] for j in indices(rm)])
        b.cx(other,lead);b.cx(cm,flag);b.push_condition(flag);b.cx(lead,dst);b.pop_condition();b.cx(cm,flag);b.cx(other,lead)
    b.free(km.array([flag]))


def phase_inverse(b,source,mask):
    """Ungated mask dot inverse(source), preserving source and classical mask."""
    info=tower(12);apply_map(b,source,info['ops'],True);apply_map(b,mask,trans_ops(12))
    a,c=source[:6],source[6:];norm=b.alloc_qubits(6);v=b.alloc_qubits(6)
    product(b,a,c,norm);linear_xor(b,a,norm,info['square']);linear_xor(b,c,norm,info['musquare'])
    def vmap():
        b.cx(c,a);classical_transpose(b,a,mask[:6],v);b.cx(c,a)
        classical_transpose(b,c,mask[6:],v)
    vmap();leaf,iv,cleanup=inverse_dot_data(b,norm,v);b.cz(leaf,iv);cleanup()
    vmap();b.free(v)
    linear_xor(b,a,norm,info['square']);linear_xor(b,c,norm,info['musquare'])
    b.del_gf2_mul(a,c,target=norm,field=km.GF2Field(6));b.free(norm)
    apply_map(b,mask,trans_ops(12),True);apply_map(b,source,info['ops'])


def apply(b,root,source,mask):
    info=tower(12);apply_map(b,source,info['ops'],True);apply_map(b,mask,trans_ops(12))
    a,c=source[:6],source[6:];norm=b.alloc_qubits(6);v=b.alloc_qubits(6)
    product(b,a,c,norm);linear_xor(b,a,norm,info['square']);linear_xor(b,c,norm,info['musquare'])
    def vmap():
        b.cx(c,a);classical_transpose(b,a,mask[:6],v);b.cx(c,a)
        classical_transpose(b,c,mask[6:],v)
    vmap();leaf,inv,cleanup=inverse_dot_data(b,norm,v)
    parity=b.alloc_qubits(1)[0];ladder=b.alloc_qubits(10);meas=b.alloc_bits(1)[0]
    for j in range(3):b.ccx(leaf[j],inv[j],parity)
    b.x(root);b.init_and(root[0],root[1],ladder[0])
    for j in range(1,10):b.init_and(ladder[j-1],root[j+1],ladder[j])
    b.ccz(ladder[-1],root[-1],parity)
    for j in reversed(range(1,10)):b.del_and(ladder[j-1],root[j+1],ladder[j])
    b.del_and(root[0],root[1],ladder[0]);b.x(root)
    b.hmr(parity,meas);b.push_condition(meas);b.cz(leaf,inv);b.pop_condition();cleanup()
    vmap();b.free(v)
    linear_xor(b,a,norm,info['square']);linear_xor(b,c,norm,info['musquare'])
    b.del_gf2_mul(a,c,target=norm,field=km.GF2Field(6));b.free(norm)
    apply_map(b,mask,trans_ops(12),True);apply_map(b,source,info['ops'])
    b.bit_store(meas,False);b.free(km.array([meas]));b.free(km.array([parity]));b.free(ladder)


def test():
    import numpy as np
    from fast_dqi import FastBuilder
    b=FastBuilder(km.CircuitBuilder());root=b.create_quantum_register(12,'root');d=b.create_quantum_register(12,'d');mask=b.create_classical_register(12,'m')
    apply(b,root,d,mask);c=b.finish_circuit();N=4096;rng=np.random.default_rng(71839)
    dd=np.arange(N);ss=rng.integers(0,4096,N);ss[::2]=0;mm=rng.integers(0,4096,N)
    bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    sim=km.Simulator(batch_size=N,seed=921);sim.use_same_registers_as(c)
    for reg,values in [(root,ss),(d,dd),(mask,mm)]:sim.write_across_shots(reg,bits(values))
    sim.do(c);gf=km.GF2Field(12)
    expected=np.array([s==0 and (gf.pow(int(d0),4094)&int(m)).bit_count()%2==1 for s,d0,m in zip(ss,dd,mm)])
    assert np.array_equal(sim.read_phase_flipped_across_shots(),expected)
    for reg,values in [(root,ss),(d,dd),(mask,mm)]:assert np.array_equal(sim.read_across_shots(reg),bits(values))
    public=set(root)|set(d);scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    assert not sim.read_across_shots(scratch).any()
    print(json.dumps(dict(status='PASS',divisors=N,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like)),flush=True)

def test_ungated():
    """Every source against every mask basis vector, plus mixed and zero masks."""
    import numpy as np
    from pathlib import Path
    from fast_dqi import FastBuilder
    b=FastBuilder(km.CircuitBuilder())
    source=b.create_quantum_register(12,'source');mask=b.create_classical_register(12,'mask')
    phase_inverse(b,source,mask);c=b.finish_circuit();gf=km.GF2Field(12)
    rng=np.random.default_rng(819291);values=np.arange(4096,dtype=np.int64)
    inverse=np.array([gf.pow(int(x),4094) for x in values])
    bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    public=set(source);scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    cases=[np.full(4096,1<<j,dtype=np.int64) for j in range(12)]
    cases += [np.zeros(4096,dtype=np.int64),rng.integers(0,4096,4096)]
    for seed,masks in enumerate(cases):
        sim=km.Simulator(batch_size=4096,seed=7291+seed);sim.use_same_registers_as(c)
        sim.write_across_shots(source,bits(values));sim.write_across_shots(mask,bits(masks));sim.do(c)
        want=np.array([int(x & y).bit_count()%2 for x,y in zip(inverse,masks)],dtype=bool)
        assert np.array_equal(sim.read_phase_flipped_across_shots(),want),seed
        assert np.array_equal(sim.read_across_shots(source),bits(values)),seed
        assert np.array_equal(sim.read_across_shots(mask),bits(masks)),seed
        assert not sim.read_across_shots(scratch).any(),seed
    out=dict(status='PASS',scope='ungated inverse parity: all field elements, all 12 mask basis vectors, zero and random masks',
             shots=4096*len(cases),Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like)
    Path('results/inverse_mask_ungated_check.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--ungated',action='store_true');args=parser.parse_args()
    test_ungated() if args.ungated else test()
