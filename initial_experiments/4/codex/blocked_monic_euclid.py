"""Small-workspace generic quotient-block Euclid update.

Copies only 2b leading coefficients to a temporary bank. Three cofactor passes
permit one full numerator product and one old-cofactor division of the tail.
Degree gaps and stopping are not handled; this is not a submission decoder.
"""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from mobius_series import BoundedArithmetic,apply as mobius
from continued_fraction_series import apply as cf
from monic_compact_euclid import step as reference_step
from compact_euclid import mul,inv,product,xor as polynomial_xor
from block_euclid_ratio import series_product,series_divide
ROOT=Path(__file__).resolve().parent


def cofactor(m,E,Bp,B,beta,alpha,pinv,inverse=False,norm=False,keep=None):
    def xor(a,b):
        for x,y in zip(a,b):E.cx(x,y)
    def mac(a,b,c):E.xor(m.prod_rt(c,a,b))
    if inverse:
        xor(B[-1],alpha);mac(alpha,B[-1],beta)
        if len(B)>1:xor(B[-2],beta)
    with m.inv_word(E,beta) as bi,m.opened(E,'block_product',12) as product:
        for i,word in enumerate(Bp):
            if keep is not None and i<len(Bp)-keep:continue
            if not inverse:m.mon_slot(E,word,beta,bi,product)
            mac(alpha,B[i],word)
            if i:xor(B[i-1],word)
            if inverse:m.mon_slot(E,word,bi,beta,product)
        if norm:m.mon_slot(E,pinv,beta if inverse else bi,bi if inverse else beta,product)
    if not inverse:
        mac(alpha,B[-1],beta)
        if len(B)>1:xor(B[-2],beta)
        xor(B[-1],alpha)


def schedule(Bp,B,Q):
    frames=[]
    for j in range(len(Q)//2):
        beta,alpha=Q[2*j:2*j+2];frames.append((Bp,B,beta,alpha));Bp,B=B,Bp+[beta,alpha]
    return frames,(Bp,B)


def apply(m,E,math,Bp,B,C,pinv,block,inverse=False,flags=None,padding=False,monomial=False,thin=False,suffix=False,linear_padding=False):
    head=list(reversed(C))[:2*block];tail=list(reversed(C)) [2*block:]
    rev=[B[0],*reversed(B)];drev=[B[0],*reversed(Bp)]
    def xor(a,b):
        for x,y in zip(a,b):E.cx(x,y)
    def head_transform(T,inv=False,record=True):
        if inv:cf(m,E,math,T,True,flags,record)
        mobius(m,E,math,rev,drev,pinv,T,inv,p_is_inverse=True)
        if not inv:cf(m,E,math,T,False,flags,record)
    with m.opened(E,'block_head',24*block) as extra:
        Q=[extra[12*j:12*j+12] for j in range(2*block)]
        frames,(finalprev,final)=schedule(Bp,B,Q);finalrev=[final[0],*reversed(final)]
        def padded_pair(index,inv,norm):
            from quantum_coprime_padding import padding_ops,apply as apply_padding
            from coprime_padding import QUADRATICS
            first,second=frames[index:index+2]
            def fake(frame,control):
                oldprev,old,be,al=frame
                for i,word in enumerate(oldprev):
                    if i:
                        for x,y in zip(old[i-1],word):E.b.ccx(control,x,y)
                if len(old)>1:
                    for x,y in zip(old[-2],be):E.b.ccx(control,x,y)
                for x,y in zip(old[-1],al):E.b.ccx(control,x,y)
            with m.opened(E,'padding_controls',2) as control:
                joint,single=control;rt=m.one_rt((flags[index],),0,(flags[index+1],),0,joint)
                E.mac(rt);E.cx(flags[index+1],single);E.cx(joint,single)
                quad=padding_ops(first[0],first[1],Q[2*index:2*index+4],True,QUADRATICS[(len(first[1])-1)//2])
                linear=padding_ops(second[0],second[1],Q[2*index+2:2*index+4],False)
                if not inv:
                    cofactor(m,E,*first,pinv,False,norm);cofactor(m,E,*second,pinv,False,norm)
                    fake(second,joint);fake(first,joint);apply_padding(E.b,quad,joint)
                    fake(second,single);apply_padding(E.b,linear,single)
                else:
                    apply_padding(E.b,linear,single,True);fake(second,single)
                    apply_padding(E.b,quad,joint,True);fake(first,joint);fake(second,joint)
                    cofactor(m,E,*second,pinv,True,norm);cofactor(m,E,*first,pinv,True,norm)
                E.cx(joint,single);E.cx(flags[index+1],single);E.free(rt)
        needed=max(len(tail)-1,2*block-1)
        def cof(inv=False,norm=False,partial=False):
            if monomial:
                from deferred_euclid_padding import fake_undo_ops,common_padding_ops
                from quantum_coprime_padding import apply as emit_ops
                assert flags is not None and not padding
                for index in (reversed(range(block)) if inv else range(block)):
                    oldprev,old,be,al=frames[index]
                    undo=fake_undo_ops(oldprev,old,be,al)
                    common=common_padding_ops(oldprev,old,[be,al],False)
                    if not inv:
                        cofactor(m,E,oldprev,old,be,al,pinv,False,norm,needed+2*index if partial else None)
                        if linear_padding:
                            from controlled_linear import apply as linear
                            linear(E.b,undo+common,flags[index])
                        else:emit_ops(E.b,undo,flags[index]);emit_ops(E.b,common,flags[index])
                    else:
                        if linear_padding:
                            from controlled_linear import apply as linear
                            linear(E.b,undo+common,flags[index],True)
                        else:emit_ops(E.b,common,flags[index],True);emit_ops(E.b,undo,flags[index],True)
                        cofactor(m,E,oldprev,old,be,al,pinv,True,norm,needed+2*index if partial else None)
            elif padding:
                assert flags is not None and block%2==0 and len(B)%2==1
                for index in (reversed(range(0,block,2)) if inv else range(0,block,2)):padded_pair(index,inv,norm)
            else:
                for index in (reversed(range(block)) if inv else range(block)):
                    cofactor(m,E,*frames[index],pinv,inv,norm,needed+2*index if partial else None)
        if not inverse:
            for source,target in zip(head,Q):xor(source,target)
            head_transform(Q)
            cof()
            if suffix:math.triangular_suffix(finalrev,1,list(reversed(C)),len(head))
            else:
                math.triangular(finalrev,1,list(reversed(C)))
                math.triangular(finalrev,1,head,True)
            cof(True,partial=thin)
            math.triangular(rev,1,tail,True)
            head_transform(head,record=False)
            for source,target in zip(Q,head):xor(source,target)
            cof(norm=True,partial=thin)
            for source,target in zip(Q,head):E.swaps(None,list(zip(source,target)))
        else:
            for source,target in zip(Q,head):E.swaps(None,list(zip(source,target)))
            cof(True,norm=True,partial=thin)
            for source,target in zip(Q,head):xor(source,target)
            head_transform(head,True,record=False)
            math.triangular(rev,1,tail)
            cof(partial=thin)
            if suffix:math.triangular_suffix(finalrev,1,list(reversed(C)),len(head),True)
            else:
                math.triangular(finalrev,1,head)
                math.triangular(finalrev,1,list(reversed(C)),True)
            cof(True)
            head_transform(Q,True)
            for source,target in zip(head,Q):xor(source,target)
    # Replace temporary-bank handles with their public destinations.
    substitution={q:t for word,target in zip(Q,head) for q,t in zip(word,target)}
    translate=lambda words:[[substitution.get(q,q) for q in word] for word in words]
    return translate(finalprev),translate(final),list(reversed(tail)),pinv


def build(N,k,block,which,armed=False):
    m=load();b=m.Builder();words=[b.qubits(12) for _ in range(N)]
    Bp=words[:k-1];B=words[k-1:2*k-1];C=words[2*k-1:-1];pinv=words[-1]
    flags=b.qubits(block) if armed else None
    sA,sB,sT=b.qubits(3);E=m.Em(b,sA,sB,sT,m.Alloc(b));math=BoundedArithmetic(m,E)
    for direction in ([False,True] if which=='both' else [which=='inverse']):apply(m,E,math,Bp,B,C,pinv,block,direction,flags)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);tag=f'blocked_monic_n{N}_k{k}_b{block}'+('_armed' if armed else '')+f'_{which}'
    tape=ROOT/f'parameters/{tag}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
    for i in range(N):kb.create_quantum_register(12,f'S{i}')
    if armed:kb.create_quantum_register(block,'zero_pivots')
    pre=ROOT/f'parameters/{tag}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=ROOT/'parameters/blocked_monic.empty';empty.write_text('');out=ROOT/f'candidates/{tag}.kmx'
    subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c


def armed_reference(k,block,bp,b,c,p):
    old_b=list(b);old_c=list(c);n=2*block;t=list(reversed(c))[:n]
    f=series_divide(series_product(t,list(reversed(b)),n),[inv(p),0]+series_product(t,list(reversed(bp)),max(0,n-2)),n)
    flags=[]
    for _ in range(block):
        zero=not f[0];flags.append(zero)
        if zero:f[0]=1
        beta=f[0];u=series_divide([beta],f,len(f));alpha=u[1];f=u[2:]
        bp,b=b,polynomial_xor(product([alpha,1],b),[mul(beta,x) for x in bp]);p=mul(p,inv(beta))
    numerator=series_product(list(reversed(b)),list(reversed(old_c)),len(old_c))
    c=list(reversed(series_divide(numerator[2*block:],list(reversed(old_b)),len(old_c)-2*block)))
    return (bp,b,c,p),flags


def test(N,k,block,shots,armed=False):
    assert 1<=k and k+block<=N//2
    rng=np.random.default_rng(318479);values=[];expected=[];flag_values=[];rejected=0
    bpidx=list(range(k-1));bidx=list(range(k-1,2*k-1));cidx=list(range(2*k-1,N-1));pidx=N-1
    Qidx=list(reversed(cidx))[:2*block];_,(newbpidx,newbidx)=schedule(bpidx,bidx,Qidx)
    while len(values)<shots:
        s=rng.integers(0,4096,N).tolist();s[-1]=int(rng.integers(1,4096));s[-2]=int(rng.integers(1,4096))
        if armed:
            if len(values)%8==0:s[2*k-1:-1]=[0]*len(cidx)
            elif len(values)%8==1:s[-2]=0
        try:
            state=(s[:k-1]+[1],s[k-1:2*k-1]+[1],s[2*k-1:-1],s[-1])
            if armed:state,flag=armed_reference(k,block,*state)
            else:
                for j in range(block):state=reference_step(k+j,*state)
        except ValueError:rejected+=1;continue
        bp,b,c,p=state;out=s.copy()
        for i,v in zip(newbpidx,bp[:-1]):out[i]=v
        for i,v in zip(newbidx,b[:-1]):out[i]=v
        for i,v in zip(cidx[:-2*block],c):out[i]=v
        out[pidx]=p;values.append(s);expected.append(out)
        if armed:flag_values.append(flag)
    values=np.array(values).T;expected=np.array(expected).T;bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    results=[];began=time.time()
    for which in ['forward','inverse','both']:
        c=build(N,k,block,which,armed);sim=km.Simulator(batch_size=shots,seed=5193);sim.use_same_registers_as(c)
        initial=expected if which=='inverse' else values;want=expected if which=='forward' else values
        for i,row in enumerate(initial):sim.write_across_shots(c.register_data[f'S{i}'],bits(row))
        if armed and which=='inverse':sim.write_across_shots(c.register_data['zero_pivots'],np.array(flag_values,dtype=bool).T)
        sim.do(c)
        for i,row in enumerate(want):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(row)),(which,i)
        if armed:assert np.array_equal(sim.read_across_shots(c.register_data['zero_pivots']),np.array(flag_values,dtype=bool).T if which=='forward' else np.zeros((block,shots),dtype=bool)),(which,'flags')
        public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        assert not sim.read_across_shots(scratch).any(),(which,'scratch')
        assert not sim.read_phase_flipped_across_shots().any(),(which,'phase')
        results.append(dict(which=which,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,extra_qubits=c.num_qubits-12*N))
    out=dict(status='PASS',scope='compact Euclid quotient block, independent reference, inverse, composition, scratch and phase; zero arming is reversible but does not decode gap/stopping paths',N=N,k=k,block=block,armed=armed,shots=shots,skipped_degree_gaps=rejected,circuits=results,seconds=time.time()-began)
    tag=f'blocked_monic_n{N}_k{k}_b{block}'+('_armed' if armed else '')
    (ROOT/f'results/{tag}_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=20);p.add_argument('--k',type=int,default=3);p.add_argument('--block',type=int,default=2);p.add_argument('--shots',type=int,default=512);p.add_argument('--armed',action='store_true');a=p.parse_args();test(a.n,a.k,a.block,a.shots,a.armed)
