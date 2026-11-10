"""Independent exact-series model and native tests for the two-pair map."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import kickmix as km
import numpy as np
import blocked_cfraction as block
from import_published import load

HERE=Path(__file__).resolve().parent
MODE='blocked_pair'
PAIRS=2
FORCE_ZERO_HEAD=False


def product(m,a,b,n):
    out=[0]*n
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<n:out[i+j]^=m.gmul(x,y)
    return out


def divide(m,a,b,n):
    out=[0]*n;inverse=m.ginv(b[0])
    for i in range(n):
        value=a[i] if i<len(a) else 0
        for j in range(1,min(i+1,len(b))):value^=m.gmul(b[j],out[i-j])
        out[i]=m.gmul(value,inverse)
    return out


def pair(m,words):
    a,u,*tail=words
    normalized=[1]+[m.gmul(x,m.ginv(u or 1)) for x in tail]
    reciprocal=divide(m,[1],normalized,len(normalized))
    b=m.gmul(u,m.ginv(a or 1));out=reciprocal[1:]
    if out:out[0]^=b
    return [a,b,*out]


def parameters(m,H):
    a,b,c,d=H;delta=b^c;epsilon=delta^d;eta=m.gmul(b,d)
    u=m.gmul(a,b);p=m.gmul(u,m.gmul(c,d));v=m.gmul(u,delta)
    w=m.gmul(u,m.gmul(delta,delta))^p
    return [1,epsilon,eta],[1,delta],p,[m.gmul(epsilon,w)^m.gmul(eta,v),m.gmul(eta,w)]


def model(m,words,inverse=False):
    H,R=list(words[:4]),list(words[4:]);n=len(R)
    if not inverse:H=pair(m,H);H=H[:2]+pair(m,H[2:])
    A,D,p,P=parameters(m,H);p=p or 1
    if n:
        if inverse:
            R=divide(m,R,A,n)
            R=product(m,D,R,n)
            R=divide(m,R,[1,*R],n)
            R=[m.gmul(p,v) for v in R]
            R=divide(m,R,D,n)
            for i,v in enumerate(P[:n]):R[i]^=v
            R=divide(m,R,A,n)
        else:
            R=product(m,A,R,n)
            for i,v in enumerate(P[:n]):R[i]^=v
            denominator=[p,*product(m,D,R,n-1)]
            R=divide(m,product(m,A,R,n),denominator,n)
    if inverse:
        def unpair(words):
            a,b,*tail=words;u=m.gmul(b,a or 1)
            if tail:tail[0]^=b
            reciprocal=divide(m,[1],[1,*tail],len(tail)+1)
            return [a,u,*[m.gmul(v,u or 1) for v in reciprocal[1:]]]
        H=H[:2]+unpair(H[2:]);H=unpair(H)
    return H+R


def build(N,inverse=False):
    m=load();m.inv_word=block.bind(m);b=m.Builder();W=[b.qubits(12) for _ in range(N)]
    E=m.Em(b,*b.qubits(3),m.Alloc(b));math=block.CubeArithmetic(m,E)
    block.block(m,E,math,W,inverse)
    tag=f'{MODE}_n{N}_'+('inverse' if inverse else 'forward')
    tape=HERE/'parameters'/f'{tag}.npy';rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);np.save(tape,rows)
    kb=km.CircuitBuilder()
    for i in range(N):kb.create_quantum_register(12,f'W{i}')
    pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=HERE/'parameters/empty';empty.write_text('');text=HERE/'candidates'/f'{tag}.kmx'
    subprocess.run([str(probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(text),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
    c=km.Circuit.from_file(str(text),format='kmx');c.to_file(str(text.with_suffix('.kmb')),format='kmb')
    return c


def probe_translator():return HERE.parent.parent/'4/codex/translate_prefix'


def check(N,shots=512):
    m=load();rng=np.random.default_rng(2026102600+N)
    values=rng.integers(0,4096,(shots,N)).tolist()
    values[:3]=[[0]*N,[1]*N,[1,0,*([0]*(N-2))]]
    values[3][:4]=[0,1,2,3]
    values[4][:4]=[1,0,2,3]
    if FORCE_ZERO_HEAD:
        for row in values[shots//2:]:row[0]=0
    expected=[model(m,row) for row in values]
    assert all(model(m,row,True)==old for row,old in zip(expected,values))
    serial_matches=0
    for row,new in zip(values,expected):
        serial=list(row)
        for i in range(PAIRS):serial=serial[:2*i]+pair(m,serial[2*i:])
        serial_matches+=serial==new
    cases=[]
    bits=lambda v:((np.asarray(v,dtype=np.uint64)[None,:]>>np.arange(12,dtype=np.uint64)[:,None])&1).astype(bool)
    for inverse in (False,True):
        c=build(N,inverse);source=expected if inverse else values;want=values if inverse else expected
        for seed in (20261026,20261027):
            sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
            for i in range(N):sim.write_across_shots(c.register_data[f'W{i}'],bits([r[i] for r in source]))
            sim.do(c)
            for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'W{i}']),bits([r[i] for r in want])),(N,inverse,i)
            assert not sim.read_phase_flipped_across_shots().any(),(N,inverse,'phase')
            assert not sim.read_across_shots(km.array([km.q(i) for i in range(12*N,c.num_qubits)])).any(),(N,inverse,'scratch')
        cases.append(dict(inverse=inverse,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=2*shots))
    report=dict(status='PASS_COMPONENT',N=N,serial_matches=serial_matches,total_series=shots,circuits=cases,
                scope='independent rational-series map; forward and inverse native, zero heads, exact phase and scratch; unsupported charts are not decoding claims')
    print(json.dumps(report),flush=True);return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,nargs='+',default=[6,14,34]);p.add_argument('--compact',action='store_true');a=p.parse_args()
    if a.compact:
        import compact_blocked_cfraction as compact
        block.block=compact.block;block.CubeArithmetic=compact.Arithmetic;MODE='compact_blocked_pair'
        def parameters(m,H):
            a,b,c,d=H;delta=b^c;eta=m.gmul(b,d);e=d or c;epsilon=delta^e
            u=m.gmul(a,b or 1);p=m.gmul(u,eta^m.gmul(delta,e));v=m.gmul(u,delta)
            w=m.gmul(u,m.gmul(delta,delta))^p
            return [1,epsilon,eta],[1,delta],p,[m.gmul(epsilon,w)^m.gmul(eta,v),m.gmul(eta,w)]
    reports=[]
    for n in a.n:
        reports.append(check(n));(HERE/f'results/{MODE}_checks.json').write_text(json.dumps(reports,indent=2)+'\n')
