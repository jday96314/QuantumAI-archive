"""Generic local Euclid transform for a compact quotient block.

Converts N syndrome words to two monic degree N/2-1,N/2 cofactors and
one inverse remainder scale. Degree gaps and stopping remain excluded.
"""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from mobius_series import BoundedArithmetic
from monic_euclid_step import apply as generic_step,inverse_keep,replace_inverse
from monic_compact_euclid import first,step
ROOT=Path(__file__).resolve().parent


def first_step(m,E,math,S,inverse=False):
    s,r=S[-2:]
    with m.opened(E,'first_inverse',12) as ri:
        inverse_keep(m,E,r,ri)
        if not inverse:
            with m.opened(E,'first_product',12) as scratch:m.mon_slot(E,s,ri,r,scratch)
        def last():
            math.product(math.spec([(1,s,2)]),math.spec([(1,ri if inverse else r,1)]),lambda tb:[S[-3][b] for b in tb])
        if inverse:last()
        for j in (reversed(range(1,len(S)-2)) if inverse else range(1,len(S)-2)):
            E.xor(m.prod_rt(S[j-1],s,S[j]))
        if not inverse:last()
        else:
            with m.opened(E,'first_product',12) as scratch:m.mon_slot(E,s,ri,r,scratch)
        replace_inverse(m,E,r,ri)


def frames(S):
    bp=[];b=[S[-2]];c=S[:-2];p=S[-1];out=[]
    for k in range(1,len(S)//2):
        out.append((bp,b,c,p));bp,b,c,p=b,bp+[p,c[-2]],c[:-2],c[-1]
    return out,(bp,b,c,p)


def build(N,which):
    m=load();b=m.Builder();words=[b.qubits(12) for _ in range(N)];sA,sB,sT=b.qubits(3)
    E=m.Em(b,sA,sB,sT,m.Alloc(b));math=BoundedArithmetic(m,E);stages,final=frames(words)
    if which in ['forward','both']:
        first_step(m,E,math,words)
        for frame in stages:generic_step(m,E,math,*frame)
    if which in ['inverse','both']:
        for frame in reversed(stages):generic_step(m,E,math,*frame,True)
        first_step(m,E,math,words,True)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);tag=f'monic_local_n{N}_{which}'
    tape=ROOT/f'parameters/{tag}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
    for i in range(N):kb.create_quantum_register(12,f'S{i}')
    pre=ROOT/f'parameters/{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=ROOT/'parameters/monic_local.empty';empty.write_text('');out=ROOT/f'candidates/{tag}.kmx'
    subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c


def test(N,shots):
    assert N>=4 and N%2==0
    rng=np.random.default_rng(193719);values=[];outputs=[];skipped=0
    stageframes,final=frames(list(range(N)));bpidx,bidx,_,pidx=final
    while len(values)<shots:
        s=rng.integers(0,4096,N).tolist()
        try:
            state=first(s)
            for k in range(1,N//2):state=step(k,*state)
        except ValueError:skipped+=1;continue
        bp,b,c,p=state;out=list(s)
        for i,v in zip(bpidx,bp[:-1]):out[i]=v
        for i,v in zip(bidx,b[:-1]):out[i]=v
        out[pidx]=p;values.append(s);outputs.append(out)
    values=np.array(values).T;outputs=np.array(outputs).T;results=[];began=time.time()
    bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    for which in ['forward','inverse','both']:
        c=build(N,which);sim=km.Simulator(batch_size=shots,seed=5312);sim.use_same_registers_as(c)
        initial=outputs if which=='inverse' else values;want=outputs if which=='forward' else values
        for i,row in enumerate(initial):sim.write_across_shots(c.register_data[f'S{i}'],bits(row))
        sim.do(c)
        for i,row in enumerate(want):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(row)),(which,i)
        public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        assert not sim.read_across_shots(scratch).any(),(which,'scratch')
        assert not sim.read_phase_flipped_across_shots().any(),(which,'phase')
        results.append(dict(which=which,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,extra_qubits=c.num_qubits-12*N))
    out=dict(status='PASS',scope='generic local Euclid quotient block, complete forward/inverse/composition, no degree gaps or decoder phase',N=N,shots=shots,skipped_degree_gaps=skipped,circuits=results,seconds=time.time()-began)
    (ROOT/f'results/monic_local_n{N}_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=16);p.add_argument('--shots',type=int,default=512);a=p.parse_args();test(a.n,a.shots)
