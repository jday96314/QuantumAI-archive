"""Quantum generic monic compressed Euclid step, with exactly N public words.

This research primitive assumes nonzero r and pinv. It does not handle degree
gaps, stopping, or decode an error word. The upstream field convention is used.
"""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from mobius_series import BoundedArithmetic
from monic_compact_euclid import step as reference_step
ROOT=Path(__file__).resolve().parent


def inverse_keep(m,E,source,out):
    wins,rts=m.chain_fwd(E,source,out)
    E.frees(rts[:-1]);m.chain_close(E,wins)


def replace_inverse(m,E,source,inverse):
    wins,rts=m.chain_open(E,inverse,source)
    E.macs(rts[:-1]);E.free(rts[-1]);E.frees(rts[:-1]);m.chain_close(E,wins)
    E.swaps(None,list(zip(source,inverse)))


def apply(m,E,math,Bp,B,C,pinv,inverse=False,boundary_flag=None):
    k=len(B);s,r=C[-2:];bulk=C[:-2]
    new=Bp+[pinv,s]
    # Leading monic coefficient is implicit. Triangular never reads A[0].
    oldrev=[B[0],*reversed(B)]
    newrev=[new[0],*reversed(new)]
    def xor(a,b):
        for x,y in zip(a,b):E.cx(x,y)
    def mac(a,b,c):E.xor(m.prod_rt(c,a,b))
    def boundary(original_r):
        if boundary_flag is not None:
            with m.opened(E,'step_boundary_mask',12) as masked:
                rt=[term for src,dst in zip(original_r,masked)
                    for term in m.one_rt((boundary_flag,),1,(src,),0,dst)]
                E.mac(rt)
                boundary_emit(masked)
                E.free(rt)
        else:boundary_emit(original_r)
    def boundary_emit(original_r):
        with m.opened(E,'step_old_s',12) as original_s:
            rt=m.prod_rt(original_s,s,original_r);E.mac(rt)
            for l,target in enumerate(bulk):
                h=len(bulk)-l
                if h<=k+1:mac(new[k+1-h],original_s,target)
                if h+1<=k+1:mac(new[k-h],original_r,target)
            E.free(rt)
    with m.opened(E,'step_rinv',12) as ri:
        inverse_keep(m,E,r,ri)
        with m.opened(E,'step_product',12) as scratch:
            if not inverse:
                m.mon_slot(E,s,ri,r,scratch)
                m.mon_slot(E,pinv,r,ri,scratch)
        if inverse:
            math.triangular(oldrev,1,list(reversed(bulk)))
            boundary(ri)
            math.triangular(newrev,1,list(reversed(bulk)),True)
            xor(B[-1],s)
            mac(s,B[-1],pinv)
            if k>1:xor(B[-2],pinv)
        else:xor(B[-1],s)
        if Bp:
            with m.inv_word(E,pinv) as beta_inv,m.opened(E,'step_product',12) as scratch:
                for i,word in enumerate(Bp):
                    if not inverse:m.mon_slot(E,word,pinv,beta_inv,scratch)
                    mac(s,B[i],word)
                    if i:xor(B[i-1],word)
                    if inverse:m.mon_slot(E,word,beta_inv,pinv,scratch)
        if not inverse:
            mac(s,B[-1],pinv)
            if k>1:xor(B[-2],pinv)
        xor(B[-1],s)
        if not inverse:
            math.triangular(newrev,1,list(reversed(bulk)))
            boundary(r)
            math.triangular(oldrev,1,list(reversed(bulk)),True)
        else:
            with m.opened(E,'step_product',12) as scratch:
                m.mon_slot(E,pinv,r,ri,scratch)
                m.mon_slot(E,s,ri,r,scratch)
        replace_inverse(m,E,r,ri)


def build(N,k,which):
    m=load();b=m.Builder();words=[b.qubits(12) for _ in range(N)]
    Bp=words[:k-1];B=words[k-1:2*k-1];C=words[2*k-1:-1];pinv=words[-1]
    sA,sB,sT=b.qubits(3);E=m.Em(b,sA,sB,sT,m.Alloc(b));math=BoundedArithmetic(m,E)
    for direction in ([False,True] if which=='both' else [which=='inverse']):
        apply(m,E,math,Bp,B,C,pinv,direction)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);tag=f'monic_step_n{N}_k{k}_{which}'
    tape=ROOT/f'parameters/{tag}.npy';np.save(tape,rows)
    kb=km.CircuitBuilder()
    for i in range(N):kb.create_quantum_register(12,f'S{i}')
    pre=ROOT/f'parameters/{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=ROOT/'parameters/monic_step.empty';empty.write_text('');out=ROOT/f'candidates/{tag}.kmx'
    subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb')
    return c


def test(N,k,shots):
    assert 1<=k<N//2
    rng=np.random.default_rng(929741);values=rng.integers(0,4096,(N,shots));values[-2]=rng.integers(1,4096,shots);values[-1]=rng.integers(1,4096,shots)
    expected=values.copy()
    for shot in range(shots):
        bp,b,c,p=reference_step(k,values[:k-1,shot].tolist()+[1],values[k-1:2*k-1,shot].tolist()+[1],values[2*k-1:-1,shot].tolist(),int(values[-1,shot]))
        expected[:k-1,shot]=b[:k-1];expected[-1,shot]=b[k-1];expected[-3,shot]=b[k]
        expected[2*k-1:-3,shot]=c;expected[-2,shot]=p
    bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    results=[];began=time.time()
    for which in ['forward','inverse','both']:
        c=build(N,k,which);sim=km.Simulator(batch_size=shots,seed=8351);sim.use_same_registers_as(c)
        initial=expected if which=='inverse' else values;want=expected if which=='forward' else values
        for i,row in enumerate(initial):sim.write_across_shots(c.register_data[f'S{i}'],bits(row))
        sim.do(c)
        for i,row in enumerate(want):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(row)),(which,i)
        public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        assert not sim.read_across_shots(scratch).any(),(which,'scratch')
        assert not sim.read_phase_flipped_across_shots().any(),(which,'phase')
        results.append(dict(which=which,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,extra_qubits=c.num_qubits-12*N))
    result=dict(status='PASS',scope='generic nonzero-pivot monic compressed Euclid step, independent classical output, inverse, composition, scratch and phase',N=N,k=k,shots=shots,circuits=results,seconds=time.time()-began)
    (ROOT/f'results/monic_step_n{N}_k{k}_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=20);p.add_argument('--k',type=int,default=3);p.add_argument('--shots',type=int,default=512);a=p.parse_args();test(a.n,a.k,a.shots)
