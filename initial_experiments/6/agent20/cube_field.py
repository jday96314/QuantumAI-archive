"""A 28-AND GF(4096) cube from seven four-AND GF(16) cubes.

Represent x=a+b*t+c*t^2, t^3=kappa over GF(16). Evaluate the degree-six
polynomial (a+b*T+c*T^2)^3 at seven GF(16) points, interpolate, and reduce
modulo T^3-kappa. All interpolation and basis changes are binary linear.
Small GF(16) cube factors are reused from round-4 cubic_trace_field.py.
"""
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import kickmix as km
from cubic_trace_field import data, cube_terms, linear


def matrix_inverse(a,g):
    n=len(a);a=[list(row)+[int(i==j) for j in range(n)] for i,row in enumerate(a)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if a[i][j])
        a[j],a[pivot]=a[pivot],a[j]
        c=g.pow(a[j][j],14);a[j]=[g.mul(v,c) for v in a[j]]
        for i in range(n):
            if i==j or not a[i][j]:continue
            c=a[i][j];a[i]=[v^g.mul(c,w) for v,w in zip(a[i],a[j])]
    return [row[n:] for row in a]


@lru_cache(None)
def host_terms():
    z=data();g=km.GF2Field(4);points=list(range(7))
    interpolation=matrix_inverse([[g.pow(x,k) for k in range(7)] for x in points],g)
    terms,lc,_=cube_terms();result={};linearcols=[0]*12
    for j,point in enumerate(points):
        evaluations=[]
        for bit in range(12):
            v=linear(z['inverse_columns'],1<<bit)
            a,b,c=(v>>(4*i)&15 for i in range(3))
            evaluations.append(a^g.mul(b,point)^g.mul(c,g.mul(point,point)))
        def output(value):
            v=0
            for degree in range(7):
                coefficient=g.mul(interpolation[degree][j],value)
                coefficient=g.mul(coefficient,g.pow(z['kappa'],degree//3))
                v^=linear(z['columns'],coefficient<<(4*(degree%3)))
            return v
        for left,right,out in terms:
            lm=sum(((left&v).bit_count()&1)<<i for i,v in enumerate(evaluations))
            rm=sum(((right&v).bit_count()&1)<<i for i,v in enumerate(evaluations))
            key=tuple(sorted((lm,rm)));result[key]=result.get(key,0)^output(out)
        for bit,v in enumerate(evaluations):linearcols[bit]^=output(linear(lc,v))
    result=[(a,b,c) for (a,b),c in result.items() if c]
    assert len(result)==28
    return tuple(result),tuple(linearcols)


def evaluate(x):
    terms,cols=host_terms();value=linear(cols,x)
    for a,b,out in terms:
        if (a&x).bit_count()%2 and (b&x).bit_count()%2:value^=out
    return value


def apply(b,source,target):
    terms,cols=host_terms();scratch=b.alloc_qubits(3);left,right,temp=scratch
    def parity(word,mask,target):
        for i,q in enumerate(word):
            if mask>>i&1:b.cx(q,target)
    for a,c,out in terms:
        parity(source,a,left);parity(source,c,right)
        b.init_and(left,right,temp)
        for j,q in enumerate(target):
            if out>>j&1:b.cx(temp,q)
        b.del_and(left,right,temp)
        parity(source,c,right);parity(source,a,left)
    for i,q in enumerate(source):
        for j,t in enumerate(target):
            if cols[i]>>j&1:b.cx(q,t)
    b.free(scratch)


if __name__=='__main__':
    import numpy as np
    g=km.GF2Field(12)
    assert all(evaluate(x)==g.pow(x,3) for x in range(4096))
    b=km.CircuitBuilder();source=b.create_quantum_register(12,'x');target=b.create_quantum_register(12,'y')
    apply(b,source,target);c=b.finish_circuit()
    xs=np.arange(4096,dtype=np.uint64);ys=(xs*119+735)%4096
    expected=ys^np.array([g.pow(int(x),3) for x in xs],dtype=np.uint64)
    bits=lambda a:((a[None,:]>>np.arange(12,dtype=np.uint64)[:,None])&1).astype(bool)
    for seed in (20261025,20261026):
        sim=km.Simulator(batch_size=4096,seed=seed);sim.use_same_registers_as(c)
        sim.write_across_shots(source,bits(xs));sim.write_across_shots(target,bits(ys));sim.do(c)
        assert np.array_equal(sim.read_across_shots(source),bits(xs))
        assert np.array_equal(sim.read_across_shots(target),bits(expected))
        assert not sim.read_phase_flipped_across_shots().any()
        assert not sim.read_across_shots(km.array([km.q(i) for i in range(24,c.num_qubits)])).any()
    path=Path('candidates/cube28.kmb');c.to_file(str(path),format='kmb')
    report=dict(status='PASS',scope='all 4096 field values, arbitrary target words, two native measurement seeds',
                shots=8192,Toffolis=c.operation_statistics().num_toffoli_like,Q=c.num_qubits,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),terms=host_terms()[0],linear_columns=host_terms()[1])
    Path('results/cube28_component.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('terms','linear_columns')}),flush=True)
