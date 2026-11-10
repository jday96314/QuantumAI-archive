"""Controlled affine cofactor padding with no arithmetic workspace.

The supplied parameter words become polynomial remainders, making each map
invertible on arbitrary states. For completed paths they hold (1,0) pairs.
This standalone prototype is not yet integrated with the compact decoder.
"""
import argparse,json
from functools import lru_cache
from pathlib import Path
import numpy as np
import kickmix as km
from fast_dqi import FastBuilder
from coprime_padding import M,OMEGA,QUADRATICS,scaled
from compact_euclid import product,xor,mul


def xor_constant(ops,c,source,target):
    for i in range(12):
        column=M.gmul(c,1<<i)
        for j in range(12):
            if column>>j&1:ops.append(('cx',source[i],target[j]))


@lru_cache(None)
def constant_map(c):
    columns=[M.gmul(c,1<<i) for i in range(12)];rows=[sum(((columns[i]>>j)&1)<<i for i in range(12)) for j in range(12)];ops=[]
    for j in range(12):
        pivot=next(i for i in range(j,12) if rows[i]>>j&1)
        if pivot!=j:rows[j],rows[pivot]=rows[pivot],rows[j];ops.append(('swap',j,pivot))
        for i in range(12):
            if i!=j and rows[i]>>j&1:rows[i]^=rows[j];ops.append(('cx',j,i))
    return tuple(reversed(ops))


def scale(ops,c,word):
    for kind,i,j in constant_map(c):ops.append((kind,word[i],word[j]))


def add_constant(ops,c,word):
    for j in range(12):
        if c>>j&1:ops.append(('x',word[j]))


def multiply(ops,words,P):
    d=len(P)-1
    for i in range(len(words)-d):
        for j in range(d):xor_constant(ops,P[j],words[i+d],words[i+j])
    for j in range(d):add_constant(ops,P[j],words[len(words)-d+j])


def permute(ops,sources,destinations):
    # Move each logical source word to the designated physical destination.
    actual={tuple(w):tuple(w) for w in sources};where={v:k for k,v in actual.items()}
    for source,destination in zip(sources,destinations):
        source,destination=tuple(source),tuple(destination);position=where[source]
        if position==destination:continue
        for a,b in zip(position,destination):ops.append(('swap',a,b))
        displaced=actual[destination];actual[position],actual[destination]=displaced,source;where[displaced]=position;where[source]=destination


def padding_ops(Bp,B,Q,quadratic=True,constant=None):
    ops=[]
    if quadratic:
        constant=QUADRATICS[0] if constant is None else constant
        P=[constant,1,1];cp=[Q[0],Q[1],*Bp];cb=[Q[2],Q[3],*B]
        add_constant(ops,1,Q[0]);add_constant(ops,1,Q[2]);multiply(ops,cp,P)
        for source,target in zip(B,cp):xor_constant(ops,1,source,target)
        add_constant(ops,1,cp[len(B)]);multiply(ops,cb,P)
        permute(ops,cp+cb,Bp+Q[:2]+B+Q[2:])
    else:
        P=[0,1];cp=[Q[1],*Bp];cb=[Q[0],*B]
        add_constant(ops,1,Q[0]);multiply(ops,cp,P)
        for word in cp:scale(ops,OMEGA,word)
        for source,target in zip(B,cp):xor_constant(ops,M.gmul(OMEGA,OMEGA),source,target)
        multiply(ops,cb,P);permute(ops,cp+cb,B+Bp+Q)
    return ops


def apply(builder,ops,flag,inverse=False):
    for op in (reversed(ops) if inverse else ops):
        if op[0]=='x':builder.cx(flag,op[1])
        elif op[0]=='cx':builder.ccx(flag,op[1],op[2])
        else:builder.cswap(flag,op[1],op[2])


def test(k,shots,quadratic=True):
    count=2*k+(3 if quadratic else 1);rng=np.random.default_rng(491379);values=rng.integers(0,4096,(count,shots));flags=rng.integers(0,2,shots,dtype=np.int64).astype(bool);expected=values.copy()
    bpidx=list(range(k-1));bidx=list(range(k-1,2*k-1));qidx=list(range(2*k-1,count))
    for shot in range(shots):
        if not flags[shot]:continue
        bp=values[bpidx,shot].tolist()+[1];b=values[bidx,shot].tolist()+[1];q=values[qidx,shot].tolist()
        if quadratic:
            P=[QUADRATICS[0],1,1];np_=xor(xor(product(P,bp),b),[q[0]^1,q[1]]);nb=xor(product(P,b),[q[2]^1,q[3]])
            pi=bpidx+qidx[:2];bi=bidx+qidx[2:]
        else:
            P=[0,1];np_=xor(scaled(xor(product(P,bp),[q[1]]),OMEGA),scaled(b,mul(OMEGA,OMEGA)));nb=xor(product(P,b),[q[0]^1]);pi=bidx;bi=bpidx+qidx
        for i,v in zip(pi,np_[:-1]):expected[i,shot]=v
        for i,v in zip(bi,nb[:-1]):expected[i,shot]=v
    bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool);results=[]
    kind='quadratic' if quadratic else 'linear';tag=f'cofactor_padding_k{k}_{kind}'
    for which in ['forward','inverse','both']:
        builder=FastBuilder(km.CircuitBuilder());words=[builder.create_quantum_register(12,f'S{i}') for i in range(count)];flag=builder.create_quantum_register(1,'flag')[0]
        ops=padding_ops(words[:k-1],words[k-1:2*k-1],words[2*k-1:],quadratic)
        if which in ['forward','both']:apply(builder,ops,flag)
        if which in ['inverse','both']:apply(builder,ops,flag,True)
        c=builder.finish_circuit();sim=km.Simulator(batch_size=shots,seed=3917);sim.use_same_registers_as(c)
        initial=expected if which=='inverse' else values;want=expected if which=='forward' else values
        for word,row in zip(words,initial):sim.write_across_shots(word,bits(row))
        sim.write_across_shots(km.array([flag]),flags[None,:]);sim.do(c)
        for i,row in enumerate(want):assert np.array_equal(sim.read_across_shots(words[i]),bits(row)),(which,i)
        assert np.array_equal(sim.read_across_shots(km.array([flag]))[0],flags),(which,'flag')
        assert not sim.read_phase_flipped_across_shots().any(),(which,'phase')
        assert c.num_qubits==12*count+1,(which,'workspace')
        c.to_file(f'candidates/{tag}_{which}.kmb',format='kmb');results.append(dict(which=which,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    out=dict(status='PASS',scope='controlled affine cofactor padding with arbitrary remainder inputs, inverse, composition and phase; not integrated with decoder',k=k,kind=kind,shots=shots,circuits=results)
    Path(f'results/{tag}_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--k',type=int,default=17);p.add_argument('--shots',type=int,default=512);p.add_argument('--linear',action='store_true');a=p.parse_args();test(a.k,a.shots,not a.linear)
