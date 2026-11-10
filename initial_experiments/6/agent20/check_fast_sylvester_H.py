"""Native comparison of optimized H preparation with the independent armed map."""
import json,time,subprocess
from pathlib import Path
import numpy as np
import kickmix as km
import fast_sylvester as fast
from sylvester_padding_screen import sylvester,prod
import check_blocked_cfraction as ref
from import_published import load
from low_workspace_inverse import bind

HERE=Path(__file__).resolve().parent
h=fast.module('agent20_native_h_reference','h_phase.py');h.sylvester=fast.apply

def build(n,both=False):
    m=load();m.inv_word=bind(m);b=m.Builder();A=[b.qubits(12) for _ in range(n)];B=[b.qubits(12) for _ in range(n+1)]
    E=m.Em(b,*b.qubits(3),m.Alloc(b));H=h.prepare(m,E,A,B)
    if both:h.prepare(m,E,A,B,True)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    tag=f'fast_sylvester_H_n{n}'+('_both' if both else '')
    tape=HERE/'parameters'/f'{tag}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
    for name,count in (('A',n),('B',n+1)):
        for i in range(count):kb.create_quantum_register(12,f'{name}{i}')
    pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=HERE/'parameters/empty';empty.write_text('');out=HERE/'candidates'/f'{tag}.kmx'
    subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c,H

def model(m,a,b):
    sq=lambda x:m.gpow(x,2048)
    C=list(map(sq,a[::2]));D=list(map(sq,a[1::2]));E=list(map(sq,b[::2]));O=list(map(sq,b[1::2]));t=len(a)//2
    out=sylvester(E,O,D,C[:-1])+[C[-1]]
    for i,x in enumerate(O[:-1]):out[t+i]^=m.gmul(C[-1],x)
    direct=[x^y for x,y in zip(prod(C,O,len(a)),prod(D,E,len(a)))]
    return out,direct

def test(n,shots=512):
    start=time.time();m=load();rng=np.random.default_rng(2026110500+n);reports=[]
    for both in (False,True):
        aa=[];bb=[];wants=[];correct=0
        for j in range(shots):
            d=j%(n+1);padding=n-d
            if both:
                a=rng.integers(0,4096,n).tolist();b=rng.integers(0,4096,n).tolist()+[1]
                want=a
            else:
                roots=rng.choice(np.arange(1,4096),d,replace=False).tolist();b=[1]
                for root in roots:b=prod(b,[root,1],len(b)+1)
                b=[0]*padding+b;a=[0]*padding+rng.integers(0,4096,d).tolist()
                want,direct=model(m,a,b);correct+=want==direct
            aa.append(a);bb.append(b);wants.append(want)
        c,H=build(n,both);bits=lambda a:((np.asarray(a,dtype=np.uint64)[None,:]>>np.arange(12,dtype=np.uint64)[:,None])&1).astype(bool)
        for seed in (20261105,20261106):
            sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
            for name,values in (('A',aa),('B',bb)):
                for i in range(len(values[0])):sim.write_across_shots(c.register_data[f'{name}{i}'],bits([r[i] for r in values]))
            sim.do(c)
            targets=[c.register_data[f'A{i}'] for i in range(n)] if both else [km.array([km.q(q) for q in word]) for word in H]
            for i,word in enumerate(targets):assert np.array_equal(sim.read_across_shots(word),bits([r[i] for r in wants])),(n,both,'A/H',i)
            for i in range(n+1):assert np.array_equal(sim.read_across_shots(c.register_data[f'B{i}']),bits([r[i] for r in bb])),(n,both,'B',i)
            assert not sim.read_phase_flipped_across_shots().any(),(n,both,'phase')
            assert not sim.read_across_shots(km.array([km.q(i) for i in range(24*n+12,c.num_qubits)])).any(),(n,both,'scratch')
        reports.append(dict(both=both,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=2*shots,actual_H_correct_fixtures=correct if not both else None))
    row=dict(status='PASS_COMPONENT',degree=n,circuits=reports,seconds=time.time()-start,
             scope='forward independent totalized Sylvester model on split locators/all padding degrees; paired arbitrary inputs; phase and scratch')
    print(json.dumps(row),flush=True);return row

if __name__=='__main__':
    rows=[]
    for n in (7,17,35):
        rows.append(test(n));(HERE/'results/fast_sylvester_H_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
