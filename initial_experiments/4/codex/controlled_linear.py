"""One-Toffoli controlled transvections for an invertible affine GF(2) map.

Reduce the matrix by rank-one involutions I+u*v^T with dot(u,v)=0.
A Clifford change of basis turns each into a single controlled CNOT. The
result uses no scratch, and at most n Toffolis for an n-bit linear map.
"""
from functools import lru_cache


def bits(value):
    while value:
        q=value&-value;yield q.bit_length()-1;value^=q


@lru_cache(maxsize=256)
def synthesize(ops):
    wires=list(dict.fromkeys(q for op in ops for q in op[1:]));index={q:i for i,q in enumerate(wires)};n=len(wires)
    rows=[1<<i for i in range(n)];bias=[0]*n
    for op in ops:
        kind=op[0];a=index[op[1]]
        if kind=='x':bias[a]^=1
        else:
            b=index[op[2]]
            if kind=='cx':rows[b]^=rows[a];bias[b]^=bias[a]
            elif kind=='swap':rows[a],rows[b]=rows[b],rows[a];bias[a],bias[b]=bias[b],bias[a]
            else:raise ValueError(op)
    factors=[]
    for i in range(n):
        col=sum(((row>>i)&1)<<j for j,row in enumerate(rows))
        if col==1<<i:continue
        u=col^(1<<i)
        if col>>i&1:v=1<<i
        else:
            j=next(j for j in bits(col) if j>i);v=(1<<i)|(1<<j)
        assert (u&v).bit_count()%2==0
        parity=0
        for j in bits(v):parity^=rows[j]
        for j in bits(u):rows[j]^=parity
        factors.append((u,v))
    assert rows==[1<<i for i in range(n)]
    return wires,tuple(factors),tuple(i for i,v in enumerate(bias) if v)


def apply(builder,ops,flag,inverse=False):
    wires,factors,bias=synthesize(tuple(ops))
    def transvection(u,v):
        pivot=next(bits(u));source=next(bits(v&~(1<<pivot)))
        spread=[wires[j] for j in bits(u&~(1<<pivot))]
        parity=[wires[j] for j in bits(v&~(1<<pivot)&~(1<<source))]
        for q in spread:builder.cx(wires[pivot],q)
        for q in parity:builder.cx(q,wires[source])
        builder.ccx(flag,wires[source],wires[pivot])
        for q in reversed(parity):builder.cx(q,wires[source])
        for q in reversed(spread):builder.cx(wires[pivot],q)
    if inverse:
        for i in bias:builder.cx(flag,wires[i])
    for u,v in (factors if inverse else reversed(factors)):transvection(u,v)
    if not inverse:
        for i in bias:builder.cx(flag,wires[i])


def test():
    import json
    from pathlib import Path
    import numpy as np
    import kickmix as km
    from fast_dqi import FastBuilder
    rng=np.random.default_rng(56178);records=[]
    for n in [2,3,7,15,31,64]:
        raw=[]
        for _ in range(6*n):
            x,y=rng.choice(n,2,replace=False);kind=rng.choice(['cx','swap','x'])
            raw.append((kind,int(x)) if kind=='x' else (kind,int(x),int(y)))
        values=rng.integers(0,2,(n,512),dtype=bool);flags=rng.integers(0,2,512,dtype=bool);expected=values.copy()
        for op in raw:
            if op[0]=='x':expected[op[1]]^=flags
            elif op[0]=='cx':expected[op[2]]^=expected[op[1]]&flags
            else:
                diff=(expected[op[1]]^expected[op[2]])&flags;expected[op[1]]^=diff;expected[op[2]]^=diff
        for which in ['forward','inverse','both']:
            b=FastBuilder(km.CircuitBuilder());x=b.create_quantum_register(n,'x');f=b.create_quantum_register(1,'f')[0];ops=[(op[0],*[x[i] for i in op[1:]]) for op in raw]
            if which in ['forward','both']:apply(b,ops,f)
            if which in ['inverse','both']:apply(b,ops,f,True)
            c=b.finish_circuit();sim=km.Simulator(batch_size=512,seed=n);sim.use_same_registers_as(c)
            sim.write_across_shots(x,expected if which=='inverse' else values);sim.write_across_shots(km.array([f]),flags[None,:]);sim.do(c)
            assert np.array_equal(sim.read_across_shots(x),expected if which=='forward' else values),(n,which)
            assert np.array_equal(sim.read_across_shots(km.array([f])),flags[None,:]);assert not sim.read_phase_flipped_across_shots().any()
            records.append(dict(n=n,which=which,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    out=dict(status='PASS',scope='random invertible affine maps, quantum controls, arbitrary data, inverse/composition and phase',shots=512,rows=records)
    Path('results/controlled_linear_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':test()
