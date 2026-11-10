"""Exact native comparison of shared synthesis to literal separate folds."""
import json
from pathlib import Path
import kickmix as km
import numpy as np
import check_blocked_cfraction as check
import shared_jacobi_synthesis as shared
from import_published import load


def model(m,W):
    def one(a,b,P,Q):
        return [a,*[m.gmul(a or 1,q) for q in Q]],[b^(Q[0] if Q else 0),*[p^m.gmul(b,q)^(Q[i+1] if i+1<len(Q) else 0) for i,(p,q) in enumerate(zip(P,Q))]]
    def fold(W):
        if len(W)==2:return [W[0]],[W[1]]
        a,b,h,e=W[:4];P,Q=fold(W[4:]) if len(W)>4 else ([],[])
        if a:
            d=e^b;c=h^m.gmul(b,d);return one(a,b,*one(c,d,P,Q))
        get=lambda i:Q[i] if i<len(Q) else 0
        return [0,b,*[m.gmul(b or 1,q) for q in Q]],[e^get(0),h^m.gmul(e,get(0))^get(1),*[p^m.gmul(h,get(i))^m.gmul(e,get(i+1))^get(i+2) for i,p in enumerate(P)]]
    P,Q=fold(list(W));return [v for pair in zip(P,Q) for v in pair]


def test(N,shots=512):
    m=load();rng=np.random.default_rng(2026102900+N);values=rng.integers(0,4096,(shots,N)).tolist()
    values[:2]=[[0]*N,[1]*N]
    for i,row in enumerate(values[2:],2):
        for k in range(0,N-2,4):
            if i%2:row[k]=0
            elif i%4==2:row[k+2]=m.gmul(row[k+1],row[k+3]^row[k+1])
    want=[model(m,row) for row in values]
    bits=lambda v:((np.asarray(v,dtype=np.uint64)[None,:]>>np.arange(12,dtype=np.uint64)[:,None])&1).astype(bool)
    report=[]
    for inverse in (False,True):
        c=check.build(N,inverse);source=want if inverse else values;expected=values if inverse else want
        for seed in (20261029,20261030):
            sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
            for i in range(N):sim.write_across_shots(c.register_data[f'W{i}'],bits([row[i] for row in source]))
            sim.do(c)
            for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'W{i}']),bits([row[i] for row in expected])),(N,inverse,i)
            assert not sim.read_phase_flipped_across_shots().any(),(N,inverse,'phase')
            assert not sim.read_across_shots(km.array([km.q(i) for i in range(12*N,c.num_qubits)])).any(),(N,inverse,'scratch')
        report.append(dict(inverse=inverse,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=2*shots))
    row=dict(status='PASS_COMPONENT',N=N,circuits=report,scope='literal separate-fold model; half gap charts, quarter zero second pivots with common factors; exact values, phase, scratch, forward and inverse')
    print(json.dumps(row),flush=True);return row


if __name__=='__main__':
    check.MODE='shared_jacobi_synthesis';check.block.block=lambda m,E,math,W,inverse=False:shared.synthesize(m,E,W,inverse)
    rows=[]
    for N in (4,6,14,34):
        rows.append(test(N));Path('results/shared_synthesis_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
