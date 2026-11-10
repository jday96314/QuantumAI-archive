"""Verify quantum completion-padded block sequence against independent algebra."""
import argparse,json,time
from pathlib import Path
import numpy as np
import kickmix as km
from blocked_euclid_sequence import frames,build
from coprime_padding import block_step
from monic_compact_euclid import first
from compact_euclid import mul


def test(N,width,shots):
    assert N//2%2==1 and width%2==0
    rng=np.random.default_rng(719413);stages,final=frames(list(range(N)),width);values=[];expected=[];flag_values=[];gaps=0
    while len(values)<shots:
        weight=len(values)%(N//2+1);positions=rng.choice(np.arange(1,4096),weight,replace=False);errors=rng.integers(1,4096,weight);s=[0]*N
        for position,error in zip(positions,errors):
            power=int(error)
            for j in range(N):s[j]^=power;power=mul(power,int(position))
        a=s.copy();zero=not a[-1]
        if zero:a[-1]=1
        state=first(a);flags=[zero]
        try:
            for k,count,*_ in stages:
                state,ff,_,_=block_step(k,count,*state,quadratic=True);flags.extend(ff)
        except ValueError:gaps+=1;continue
        bp,b,c,p=state;out=s.copy();bpidx,bidx,_,pidx=final
        for i,v in zip(bpidx,bp[:-1]):out[i]=v
        for i,v in zip(bidx,b[:-1]):out[i]=v
        out[pidx]=p;values.append(s);expected.append(out);flag_values.append(flags)
    values=np.array(values).T;expected=np.array(expected).T;flag_values=np.array(flag_values,dtype=bool).T;results=[];began=time.time()
    bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    for which in ['forward','inverse','both']:
        c=build(N,width,which,padding=True);sim=km.Simulator(batch_size=shots,seed=6791);sim.use_same_registers_as(c)
        initial=expected if which=='inverse' else values;want=expected if which=='forward' else values
        if which=='both':
            # Independent random input words and forced zero leaders also have
            # to round-trip, including paths outside the decoder's gap model.
            initial=rng.integers(0,4096,(N,shots));initial[-1,::3]=0;initial[:,::17]=0;want=initial.copy()
        for i,row in enumerate(initial):sim.write_across_shots(c.register_data[f'S{i}'],bits(row))
        if which=='inverse':sim.write_across_shots(c.register_data['zero_pivots'],flag_values)
        sim.do(c)
        for i,row in enumerate(want):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(row)),(which,i)
        assert np.array_equal(sim.read_across_shots(c.register_data['zero_pivots']),flag_values if which=='forward' else np.zeros_like(flag_values)),(which,'flags')
        public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        assert not sim.read_across_shots(scratch).any(),(which,'scratch')
        assert not sim.read_phase_flipped_across_shots().any(),(which,'phase')
        results.append(dict(which=which,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    out=dict(status='PASS',scope='completion-padded full quantum compact-state sequence versus independent algebra on all weights; arbitrary gap-input composition restores inputs, flags, scratch and phase; no Wronskian or decoder phase',N=N,block=width,shots=shots,skipped_degree_gaps=gaps,circuits=results,seconds=time.time()-began)
    Path(f'results/padded_sequence_n{N}_b{width}_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=22);p.add_argument('--block',type=int,default=4);p.add_argument('--shots',type=int,default=512);a=p.parse_args();test(a.n,a.block,a.shots)
