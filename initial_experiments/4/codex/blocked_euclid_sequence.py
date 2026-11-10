"""Compose armed compact Euclid blocks through a full fixed-length transform.

This verifies state values and reversibility, including completion and gap
inputs. It does not implement unpadding, the Wronskian map, or decoder phase.
"""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from mobius_series import BoundedArithmetic
from monic_euclid_local import first_step
from monic_compact_euclid import first
from blocked_monic_euclid import apply as block_apply,schedule,armed_reference
from compact_euclid import mul
ROOT=Path(__file__).resolve().parent


def frames(S,width):
    bp=[];b=[S[-2]];c=S[:-2];p=S[-1];out=[];k=1
    while k<len(S)//2:
        count=min(width,len(S)//2-k);out.append((k,count,bp,b,c,p))
        q=list(reversed(c))[:2*count];_,(bp,b)=schedule(bp,b,q);c=c[:-2*count];k+=count
    return out,(bp,b,c,p)


def build(N,width,which,padding=False):
    m=load();b=m.Builder();words=[b.qubits(12) for _ in range(N)];flags=b.qubits(N//2);sA,sB,sT=b.qubits(3)
    E=m.Em(b,sA,sB,sT,m.Alloc(b));math=BoundedArithmetic(m,E);stages,_=frames(words,width)
    if which in ['forward','both']:
        with m.eq_test(E,words[-1],0) as top:E.cx(top,flags[0])
        E.cx(flags[0],words[-1][0]);first_step(m,E,math,words)
        for k,count,bp,bb,c,p in stages:block_apply(m,E,math,bp,bb,c,p,count,flags=flags[k:k+count],padding=padding)
    if which in ['inverse','both']:
        for k,count,bp,bb,c,p in reversed(stages):block_apply(m,E,math,bp,bb,c,p,count,True,flags[k:k+count],padding=padding)
        first_step(m,E,math,words,True);E.cx(flags[0],words[-1][0])
        with m.eq_test(E,words[-1],0) as top:E.cx(top,flags[0])
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);tag=f'blocked_sequence_n{N}_b{width}'+('_padding' if padding else '')+f'_{which}'
    tape=ROOT/f'parameters/{tag}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
    for i in range(N):kb.create_quantum_register(12,f'S{i}')
    kb.create_quantum_register(N//2,'zero_pivots');pre=ROOT/f'parameters/{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=ROOT/'parameters/blocked_sequence.empty';empty.write_text('');out=ROOT/f'candidates/{tag}.kmx'
    subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c


def test(N,width,shots):
    rng=np.random.default_rng(619714);values=[];expected=[];flag_values=[];stages,final=frames(list(range(N)),width)
    for shot in range(shots):
        if shot%2:
            s=rng.integers(0,4096,N).tolist()
            if shot%16==1:s[-1]=0
        else:
            weight=(shot//2)%(N//2+1);positions=rng.choice(np.arange(1,4096),weight,replace=False);errors=rng.integers(1,4096,weight);s=[0]*N
            for position,error in zip(positions,errors):
                power=int(error)
                for j in range(N):s[j]^=power;power=mul(power,int(position))
        armed=s.copy();zero=not armed[-1]
        if zero:armed[-1]=1
        state=first(armed);flags=[zero]
        for k,count,*_ in stages:
            state,newflags=armed_reference(k,count,*state);flags.extend(newflags)
        bp,b,c,p=state;out=s.copy();bpidx,bidx,_,pidx=final
        for i,v in zip(bpidx,bp[:-1]):out[i]=v
        for i,v in zip(bidx,b[:-1]):out[i]=v
        out[pidx]=p;values.append(s);expected.append(out);flag_values.append(flags)
    values=np.array(values).T;expected=np.array(expected).T;flag_values=np.array(flag_values,dtype=bool).T;results=[];began=time.time()
    bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    for which in ['forward','inverse','both']:
        c=build(N,width,which);sim=km.Simulator(batch_size=shots,seed=7164);sim.use_same_registers_as(c)
        initial=expected if which=='inverse' else values;want=expected if which=='forward' else values
        for i,row in enumerate(initial):sim.write_across_shots(c.register_data[f'S{i}'],bits(row))
        if which=='inverse':sim.write_across_shots(c.register_data['zero_pivots'],flag_values)
        sim.do(c)
        for i,row in enumerate(want):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(row)),(which,i)
        assert np.array_equal(sim.read_across_shots(c.register_data['zero_pivots']),flag_values if which=='forward' else np.zeros_like(flag_values)),(which,'flags')
        public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        assert not sim.read_across_shots(scratch).any(),(which,'scratch')
        assert not sim.read_phase_flipped_across_shots().any(),(which,'phase')
        results.append(dict(which=which,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    out=dict(status='PASS',scope='full armed compact-state transform and inverse, including zero syndromes, all weights and gap inputs; no unpadding, Wronskian or decoder phase',N=N,block=width,shots=shots,zero_pivot_events=int(flag_values.sum()),circuits=results,seconds=time.time()-began)
    (ROOT/f'results/blocked_sequence_n{N}_b{width}_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=20);p.add_argument('--block',type=int,default=4);p.add_argument('--shots',type=int,default=512);a=p.parse_args();test(a.n,a.block,a.shots)
