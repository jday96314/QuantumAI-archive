"""Native sampled reversibility tests; these do not accept a full decoder."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
import kickmix as km
from import_published import load
from components.trimmed_middle import TrimmedArithmetic
from components.budget_inverse import bind
from blocked_monic_euclid import apply as block,schedule
from monomial_serial_euclid import apply as serial
from check_deferred_padding import monomial_reference

HERE=Path(__file__).resolve().parent
R4=HERE.parents[1]/'4/codex'

def native(m,b,registers,tag):
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    base=HERE/'candidates'/tag;tape=base.with_suffix('.npy');np.save(tape,rows)
    kb=km.CircuitBuilder()
    for name,size in registers:kb.create_quantum_register(size,name)
    pre=base.with_suffix('.pre.kmx');kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=base.with_suffix('.empty');empty.write_text('');out=base.with_suffix('.kmx')
    subprocess.run([str(R4/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,stdout=subprocess.DEVNULL)
    c=km.Circuit.from_file(str(out),format='kmx');dest=base.with_suffix('.kmb');c.to_file(str(dest),format='kmb')
    # KMB is the reviewable circuit artifact. Remove redundant large tapes.
    tape.unlink();out.unlink();pre.unlink();empty.unlink()
    return c,dict(Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,
                  peak_live=b.peak_qubits,sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),artifact=str(dest))

def bits(v):return ((v[None,:]>>np.arange(12)[:,None])&1).astype(bool)

def simulate(c,initial,expected,seed):
    shots=next(iter(initial.values())).shape[1]
    sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
    for name,value in initial.items():sim.write_across_shots(c.register_data[name],value)
    sim.do(c)
    for name,value in expected.items():assert np.array_equal(sim.read_across_shots(c.register_data[name]),value),name
    public={q for reg in c.register_data.values() for q in reg}
    scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    assert not sim.read_across_shots(scratch).any(),'scratch'
    assert not sim.read_phase_flipped_across_shots().any(),'phase'

def middle(shots):
    rng=np.random.default_rng(725938);records=[]
    # Includes selected early/late tail dimensions, Frobenius q>1, and
    # suffix offsets that require trimmed transposed convolution.
    for degree,length,q,start in [(3,68,1,6),(22,32,1,6),(7,24,2,11),(4,21,4,9)]:
        m=load();a=rng.integers(0,4096,(degree,shots));t=rng.integers(0,4096,(length,shots))
        a[:,0]=0;t[:,1]=0;want=t.copy()
        for i in range(start,length):
            for j in range(i):
                h=i-j
                if h%q==0 and 1<=h//q<=degree:
                    want[i]^=np.array([m.gmul(m.gpow(int(x),q),int(y)) for x,y in zip(a[h//q-1],t[j])])
        for mode in ['forward','inverse','both','inverse_forward']:
            b=m.Builder();A=[b.qubits(12) for _ in range(degree)];T=[b.qubits(12) for _ in range(length)]
            E=m.Em(b,*b.qubits(3),m.Alloc(b));math=TrimmedArithmetic(m,E)
            if mode in ['inverse','inverse_forward']:math.triangular_suffix([A[0],*A],q,T,start,True)
            if mode in ['forward','both','inverse_forward']:math.triangular_suffix([A[0],*A],q,T,start)
            if mode=='both':math.triangular_suffix([A[0],*A],q,T,start,True)
            tag=f'middle_d{degree}_n{length}_q{q}_s{start}_{mode}'
            c,record=native(m,b,[(f'A{i}',12) for i in range(degree)]+[(f'T{i}',12) for i in range(length)],tag)
            inputs={f'A{i}':bits(row) for i,row in enumerate(a)}
            inputs.update({f'T{i}':bits(row) for i,row in enumerate(want if mode=='inverse' else t)})
            expected={f'A{i}':bits(row) for i,row in enumerate(a)}
            expected.update({f'T{i}':bits(row) for i,row in enumerate(want if mode=='forward' else t)})
            simulate(c,inputs,expected,61731);record.update(degree=degree,length=length,q=q,start=start,mode=mode);records.append(record)
            print(json.dumps(record),flush=True)
    return records

def blocks(shots):
    rng=np.random.default_rng(819482);records=[]
    for k,count,is_serial in [(1,3,False),(19,3,False),(22,1,True),(34,1,True)]:
        n=70;vals=[];wants=[];flags=[]
        m=load()
        while len(vals)<shots:
            s=rng.integers(0,4096,n).tolist();s[-1]=int(rng.integers(1,4096))
            if len(vals)%8==0:s[2*k-1:-1]=[0]*(n-2*k)
            try:state,ff=monomial_reference(k,count,s[:k-1]+[1],s[k-1:2*k-1]+[1],s[2*k-1:-1],s[-1])
            except ValueError:continue
            vals.append(s);flags.append(ff)
            bp,bb,c,p=state;out=s.copy();inds=list(range(n));Q=list(reversed(inds[2*k-1:-1]))[:2*count]
            _,(ibp,ibb)=schedule(inds[:k-1],inds[k-1:2*k-1],Q)
            for i,v in zip(ibp,bp[:-1]):out[i]=v
            for i,v in zip(ibb,bb[:-1]):out[i]=v
            for i,v in zip(inds[2*k-1:-1-2*count],c):out[i]=v
            out[-1]=p;wants.append(out)
        values=np.array(vals).T;expected=np.array(wants).T;flagvalues=np.array(flags,dtype=bool).T
        for mode in ['forward','inverse','both']:
            m=load();b=m.Builder();W=[b.qubits(12) for _ in range(n)];F=b.qubits(k+count)
            E=m.Em(b,*b.qubits(3),m.Alloc(b));m.inv_word=bind(m);math=TrimmedArithmetic(m,E);math.virtual_zero=True
            bp,bb,c,p=W[:k-1],W[k-1:2*k-1],W[2*k-1:-1],W[-1]
            for inv in ([False,True] if mode=='both' else [mode=='inverse']):
                if is_serial:serial(m,E,math,bp,bb,c,p,F[k],inv,linear_padding=True)
                else:block(m,E,math,bp,bb,c,p,count,inv,flags=F[k:],monomial=True,suffix=True,linear_padding=True)
            tag=f'block_k{k}_b{count}_{mode}'
            circuit,record=native(m,b,[(f'S{i}',12) for i in range(n)]+[('flags',k+count)],tag)
            inputs={f'S{i}':bits(row) for i,row in enumerate(expected if mode=='inverse' else values)}
            outputs={f'S{i}':bits(row) for i,row in enumerate(expected if mode=='forward' else values)}
            ff=np.vstack([np.zeros((k,shots),dtype=bool),flagvalues])
            inputs['flags']=ff if mode=='inverse' else np.zeros_like(ff)
            outputs['flags']=ff if mode=='forward' else np.zeros_like(ff)
            simulate(circuit,inputs,outputs,97138);assert record['peak_live']<=961
            record.update(k=k,block=count,mode=mode,forced_zero_tail_cases=(shots+7)//8);records.append(record)
            print(json.dumps(record),flush=True)
    return records

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('kind',choices=['middle','blocks']);ap.add_argument('--shots',type=int,default=512);args=ap.parse_args()
    start=time.time();rows=(middle if args.kind=='middle' else blocks)(args.shots)
    report=dict(status='PASS',scope='native component tests only; independent polynomial output, inverse and composition, preserved inputs, clean scratch and zero phase',shots=args.shots,rows=rows,seconds=time.time()-start)
    (HERE/f'results/native_{args.kind}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'}),flush=True)
