"""Expose variable-degree cofactors in fixed banks by rotating common x factors.

Completion flags encode a unary suffix. Their binary count is a linear XOR
map on that promise. Keeping the leading monic coefficients explicit permits
both banks to be rotated down without any irreducible padding or H map.
"""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km


def degree_count(b,flags,out):
    n=len(flags)
    for k,flag in enumerate(flags):
        mask=(n-k)^(n-k-1)
        for i,q in enumerate(out):
            if mask>>i&1:b.cx(flag,q)


def rotate(b,words,control,amount,inverse=False):
    n=len(words);amount=(-amount if inverse else amount)%n
    # Desired output[i]=old[(i+amount)%n]. A cycle is emitted as swaps from
    # its first index in reverse order, independent of coefficient values.
    seen=set()
    for start in range(n):
        if start in seen:continue
        cycle=[];i=start
        while i not in seen:seen.add(i);cycle.append(i);i=(i+amount)%n
        for i in reversed(cycle[1:]):b.cswap(control,words[start],words[i])


def apply(b,Bp,B,flags,counter,inverse=False):
    if not inverse:degree_count(b,flags,counter)
    for i in (reversed(range(len(counter))) if inverse else range(len(counter))):
        for bank in [Bp,B]:rotate(b,bank,counter[i],1<<i,inverse)
    if inverse:degree_count(b,flags,counter)


def test(n=35,shots=512):
    b=km.CircuitBuilder();bp=[b.create_quantum_register(12,f'P{i}') for i in range(n)];bb=[b.create_quantum_register(12,f'B{i}') for i in range(n+1)];flags=b.create_quantum_register(n,'flags');counter=b.alloc_qubits(n.bit_length())
    apply(b,bp,bb,flags,counter);forward=b.finish_circuit()
    b=km.CircuitBuilder();bp2=[b.create_quantum_register(12,f'P{i}') for i in range(n)];bb2=[b.create_quantum_register(12,f'B{i}') for i in range(n+1)];ff=b.create_quantum_register(n,'flags');cc=b.alloc_qubits(n.bit_length());apply(b,bp2,bb2,ff,cc);apply(b,bp2,bb2,ff,cc,True);both=b.finish_circuit()
    rng=np.random.default_rng(92801);pv=rng.integers(0,4096,(n,shots));bv=rng.integers(0,4096,(n+1,shots));weight=np.arange(shots)%(n+1);fv=np.arange(n)[:,None]>=weight
    bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    for which,c in [('forward',forward),('both',both)]:
        sim=km.Simulator(batch_size=shots,seed=1);sim.use_same_registers_as(c)
        for prefix,values in [('P',pv),('B',bv)]:
            for j,v in enumerate(values):sim.write_across_shots(c.register_data[f'{prefix}{j}'],bits(v))
        sim.write_across_shots(c.register_data['flags'],fv);sim.do(c)
        for prefix,values in [('P',pv),('B',bv)]:
            if which=='forward':want=np.stack([np.roll(values[:,s],-(n-int(weight[s]))%len(values)) for s in range(shots)],1)
            else:want=values
            for j,v in enumerate(want):assert np.array_equal(sim.read_across_shots(c.register_data[f'{prefix}{j}']),bits(v)),(which,prefix,j)
        assert np.array_equal(sim.read_across_shots(c.register_data['flags']),fv)
        count=sum((sim.read_across_shots(counter).astype(int)<<np.arange(len(counter))[:,None]))
        assert np.array_equal(count,n-weight if which=='forward' else np.zeros(shots,dtype=int))
        assert not sim.read_phase_flipped_across_shots().any()
    out=dict(status='PASS',scope='cofactor coefficient rotations with unary completion flags; all degrees, arbitrary polynomial values, composition and phase',degree=n,shots=shots,Q=forward.num_qubits,T=forward.operation_statistics().num_toffoli_like)
    Path('results/compact_monomial_unpad_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':test()
