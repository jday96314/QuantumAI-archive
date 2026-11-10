"""Native H repair tests with forced first gaps and every padding degree."""
import json,time,subprocess
from pathlib import Path
import numpy as np
import kickmix as km
import repaired_sylvester as repair_module
from sylvester_padding_screen import sylvester,prod
import check_blocked_cfraction as ref
from import_published import load
from low_workspace_inverse import bind

HERE=Path(__file__).resolve().parent
h=repair_module.fast.module('agent20_repaired_h_test','h_phase.py')

def build(n,both=False):
    m=load();m.inv_word=bind(m);b=m.Builder();A=[b.qubits(12) for _ in range(n)];B=[b.qubits(12) for _ in range(n+1)];flags=b.qubits(n)
    E=m.Em(b,*b.qubits(3),m.Alloc(b));repair=repair_module.Repair(m,E,flags);h.sylvester=repair.apply
    H=h.prepare(m,E,A,B)
    if both:h.prepare(m,E,A,B,True);repair.close()
    held=set(repair.tags+repair.gamma+repair.count) if not both else set()
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    tag=f'repaired_sylvester_H_n{n}'+('_both' if both else '')
    tape=HERE/'parameters'/f'{tag}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
    for name,count in (('A',n),('B',n+1)):
        for i in range(count):kb.create_quantum_register(12,f'{name}{i}')
    kb.create_quantum_register(n,'flags')
    pre=HERE/'parameters'/f'{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=HERE/'parameters/empty';empty.write_text('');out=HERE/'candidates'/f'{tag}.kmx'
    subprocess.run([str(ref.probe_translator()),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c,H,held,repair.tags

def model(m,a,b,padding):
    sq=lambda x:m.gpow(x,2048);C=list(map(sq,a[::2]));D=list(map(sq,a[1::2]));E=list(map(sq,b[::2]));O=list(map(sq,b[1::2]));t=len(a)//2
    out=sylvester(E,O,D,C[:-1],[False,padding//2])+[C[-1]]
    for i,x in enumerate(O[:-1]):out[t+i]^=m.gmul(C[-1],x)
    direct=[x^y for x,y in zip(prod(C,O,len(a)),prod(D,E,len(a)))]
    return out,direct

def test(n,shots=512):
    start=time.time();m=load();rng=np.random.default_rng(2026110600+n);reports=[];t=n//2
    for both in (False,True):
        aa=[];bb=[];flagvals=[];wants=[];forced=0;rejected=0
        for j in range(shots):
            d=j%(n+1);padding=n-d;k=padding//2
            if both:
                a=rng.integers(0,4096,n).tolist();b=rng.integers(0,4096,n).tolist()+[1]
                flags=rng.integers(0,2,n).astype(bool).tolist();want=a
            else:
                force=bool(j%2==0 and k<=t-2)
                while True:
                    if force:
                        O=[0]*k+rng.integers(1,4096,t-k).tolist()+[1]
                        R=[0]*k+rng.integers(1,4096,t-k-1).tolist()+[0]
                        q=int(rng.integers(1,4096))
                        if padding%2:R[k]=m.gmul(q,O[k])
                        E=[m.gmul(q,x)^(R[i] if i<t else 0) for i,x in enumerate(O)]
                        b=[v for e,o in zip(E,O) for v in (m.gmul(e,e),m.gmul(o,o))]
                    else:
                        roots=rng.choice(np.arange(1,4096),d,replace=False).tolist();b=[1]
                        for root in roots:b=prod(b,[root,1],len(b)+1)
                        b=[0]*padding+b
                    a=[0]*padding+rng.integers(0,4096,d).tolist()
                    candidate,direct=model(m,a,b,padding)
                    if candidate==direct:want=direct;break
                    rejected+=1
                forced+=force;flags=[i>=d for i in range(n)]
            aa.append(a);bb.append(b);flagvals.append(flags);wants.append(want)
        c,H,held,tags=build(n,both);bits=lambda a:((np.asarray(a,dtype=np.uint64)[None,:]>>np.arange(12,dtype=np.uint64)[:,None])&1).astype(bool)
        for seed in (20261106,20261107):
            sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
            for name,values in (('A',aa),('B',bb)):
                for i in range(len(values[0])):sim.write_across_shots(c.register_data[f'{name}{i}'],bits([r[i] for r in values]))
            sim.write_across_shots(c.register_data['flags'],np.asarray(flagvals,dtype=bool).T);sim.do(c)
            targets=[c.register_data[f'A{i}'] for i in range(n)] if both else [km.array([km.q(q) for q in word]) for word in H]
            for i,word in enumerate(targets):assert np.array_equal(sim.read_across_shots(word),bits([r[i] for r in wants])),(n,both,'A/H',i)
            for i in range(n+1):assert np.array_equal(sim.read_across_shots(c.register_data[f'B{i}']),bits([r[i] for r in bb])),(n,both,'B',i)
            assert np.array_equal(sim.read_across_shots(c.register_data['flags']),np.asarray(flagvals,dtype=bool).T)
            assert not sim.read_phase_flipped_across_shots().any(),(n,both,'phase')
            scratch=km.array([km.q(i) for i in range(25*n+12,c.num_qubits) if i not in held])
            assert not sim.read_across_shots(scratch).any(),(n,both,'scratch')
            if not both:
                tagbits=sim.read_across_shots(km.array([km.q(i) for i in tags])).astype(int)
                assert np.max((tagbits*(1<<np.arange(len(tags)))[:,None]).sum(0))<=n//2
        reports.append(dict(both=both,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=2*shots,
                            forced_first_gaps=forced,fixtures_outside_single_gap_domain_rejected=rejected))
    row=dict(status='PASS_COMPONENT',degree=n,circuits=reports,seconds=time.time()-start,
             scope='true H including forced first gaps/all padding degrees, intentional gamma/tag garbage retained forward; arbitrary paired inverse incl non-unary flags, phase and clean scratch')
    print(json.dumps(row),flush=True);return row

if __name__=='__main__':
    rows=[]
    for n in (7,17,35):
        rows.append(test(n));(HERE/'results/repaired_sylvester_H_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
