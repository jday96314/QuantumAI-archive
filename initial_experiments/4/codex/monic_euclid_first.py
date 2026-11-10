"""First compact, monic Euclidean step in the syndrome array itself.

For S of degree N-1 with leading a, v=S[N-2]/a, the high remainder of
a*x^N+(x+v)*S is S[j-1]+v*S[j]. Retain coefficients 1..N-2, v, and 1/a.
The omitted constant is implied by divisibility by x+v. A single flag arms
a zero leading coefficient reversibly; all input states have an exact inverse.
This is a research primitive, not a complete decoder.
"""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from fast_dqi import FastBuilder
from tower_inverse import tower,apply_map,linear_xor,del_inverse
from stream_forney import word_mul,K,GF
from inverse_mask_kernel import phase_inverse

def inverse_keep(b,source):
    out,scaffold=b.init_gf2_inverse_with_scaffold(source,field=GF)
    # Erase the norm-chain scaffold while retaining the computed inverse.
    norm,iv,inner=scaffold[:6],scaffold[6:12],scaffold[12:]
    del_inverse(b,norm,iv,inner)
    info=tower(12);apply_map(b,source,info['ops'],True);a,c=source[:6],source[6:]
    linear_xor(b,a,norm,info['square']);linear_xor(b,c,norm,info['musquare'])
    b.del_gf2_mul(a,c,target=norm,field=km.GF2Field(6));apply_map(b,source,info['ops']);b.free(scaffold)
    return out

def swap_inverse(b,original,inverse):
    mask=b.alloc_bits(12);b.hmr(original,mask);phase_inverse(b,inverse,mask)
    b.bit_store(mask,False);b.free(mask);b.swap(original,inverse);b.free(inverse)

def mac(b,a,c,target,square=False):
    if square:
        from tower_inverse import map_ops
        ops=map_ops(tuple(GF.mul(1<<i,1<<i) for i in range(12)));apply_map(b,a,ops)
    K._mul(b,a,c,target)
    if square:apply_map(b,a,ops,True)

def apply(b,S,armed,inverse=False):
    n=len(S);assert n>=4
    if not inverse:
        b.flip_if_equal(S[-1],0,target=armed);b.cx(armed,S[-1][0])
    reciprocal=inverse_keep(b,S[-1])
    b.info_tag('euclid_inv_ready')
    if not inverse:
        word_mul(b,S[-2],reciprocal,S[-1])
        for j in range(1,n-2):mac(b,S[-2],S[j],S[j-1])
        mac(b,S[-2],S[-1],S[-3],True)
    else:
        mac(b,S[-2],reciprocal,S[-3],True)
        for j in reversed(range(1,n-2)):mac(b,S[-2],S[j],S[j-1])
        b.info_tag('euclid_remainder_restored')
        word_mul(b,S[-2],reciprocal,S[-1])
        b.info_tag('euclid_sub_restored')
    swap_inverse(b,S[-1],reciprocal)
    b.info_tag('euclid_lead_restored')
    if inverse:
        b.cx(armed,S[-1][0]);b.flip_if_equal(S[-1],0,target=armed)

def test(n=70,shots=512):
    rng=np.random.default_rng(678121);values=rng.integers(0,4096,(n,shots),dtype=np.int64);values[-1,::8]=0
    def build(which):
        b=FastBuilder(km.CircuitBuilder());S=[b.create_quantum_register(12,f'S{i}') for i in range(n)];a=b.create_quantum_register(1,'armed')[0]
        if which in ['forward','both']:apply(b,S,a)
        if which in ['inverse','both']:apply(b,S,a,True)
        return b.finish_circuit(),S,a
    c,S,armed=build('forward');ci,_,_=build('inverse');sim=km.Simulator(batch_size=shots,seed=223);sim.use_same_registers_as(c)
    bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    for reg,row in zip(S,values):sim.write_across_shots(reg,bits(row))
    sim.do(c);lead=values[-1].copy();lead[lead==0]=1;iv=np.array([GF.pow(int(x),4094) for x in lead]);v=np.array([GF.mul(int(x),int(y)) for x,y in zip(values[-2],iv)])
    want=values.copy();want[-1]=iv;want[-2]=v
    for j in range(1,n-1):want[j-1]=values[j-1]^np.array([GF.mul(int(x),int(y)) for x,y in zip(v,values[j])])
    for i,row in enumerate(want):assert np.array_equal(sim.read_across_shots(S[i]),bits(row)),('forward',i)
    assert np.array_equal(sim.read_across_shots(km.array([armed]))[0],values[-1]==0);assert not sim.read_phase_flipped_across_shots().any()
    public={q for reg in S for q in reg}|{armed};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);assert not sim.read_across_shots(scratch).any()
    # Independently compiled circuits assume fresh classical scratch IDs.
    # Test inverse on fresh scratch, and composition within one builder.
    sim=km.Simulator(batch_size=shots,seed=224);sim.use_same_registers_as(ci)
    for reg,row in zip(S,want):sim.write_across_shots(reg,bits(row))
    sim.write_across_shots(km.array([armed]),np.array([values[-1]==0]));sim.do(ci)
    for i,row in enumerate(values):assert np.array_equal(sim.read_across_shots(S[i]),bits(row)),('inverse',i)
    assert not sim.read_across_shots(km.array([armed])).any();assert not sim.read_phase_flipped_across_shots().any();assert not sim.read_across_shots(scratch).any()
    combined,_,_=build('both');sim=km.Simulator(batch_size=shots,seed=225);sim.use_same_registers_as(combined)
    for reg,row in zip(S,values):sim.write_across_shots(reg,bits(row))
    sim.do(combined)
    for i,row in enumerate(values):assert np.array_equal(sim.read_across_shots(S[i]),bits(row)),('combined',i)
    public_S={q for reg in S for q in reg};junk=km.array([km.q(i) for i in range(combined.num_qubits) if km.q(i) not in public_S])
    assert not sim.read_phase_flipped_across_shots().any();assert not sim.read_across_shots(junk).any()
    out=dict(status='PASS',scope='first monic Euclidean step and its inverse; includes armed zero pivots',N=n,shots=shots,Q=c.num_qubits,forward_T=c.operation_statistics().num_toffoli_like,inverse_T=ci.operation_statistics().num_toffoli_like)
    c.to_file(f'candidates/monic_first_{n}.kmb',format='kmb');ci.to_file(f'candidates/monic_first_{n}_inverse.kmb',format='kmb')
    Path(f'results/monic_first_{n}_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--n',type=int,default=70);a.add_argument('--shots',type=int,default=512);args=a.parse_args();test(args.n,args.shots)
