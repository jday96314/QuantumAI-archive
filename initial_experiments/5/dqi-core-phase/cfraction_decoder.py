"""Reversible C-fraction decoder in the syndrome array itself.

One pair peel uses a normalized reciprocal series. Pair synthesis overwrites
the C-fraction array with the numerator/denominator, without separate banks.
Zero pivots are armed: generic charts decode, every input still restores.
This is a research decoder; it is not an accepted DQI candidate.
"""
from contextlib import contextmanager
from pathlib import Path
import argparse,hashlib,json,subprocess
import numpy as np
import kickmix as km
from import_published import load
from low_workspace_inverse import bind
from mobius_series import BoundedArithmetic
from published_dqi import isomorphism
from fast_dqi import FastBuilder
from recursive_forney import RecursiveForney
from tower_inverse import map_ops,apply_map,linear_xor
from algebra import Field,trim
from check_algebra import peel as exact_peel

HERE=Path(__file__).resolve().parent
TRANSLATOR=HERE/'../../4/codex/translate_prefix'


@contextmanager
def armed(m,E,word):
    with m.eq_test(E,word,0) as zero:
        E.cx(zero,word[0]);yield;E.cx(zero,word[0])


def scale(m,E,words,pivot,inverse=False):
    if not words:return
    with armed(m,E,pivot),m.inv_word(E,pivot) as pi,m.opened(E,'cf_scale_product',12) as t:
        for word in words:m.mon_slot(E,word,pi if inverse else pivot,pivot if inverse else pi,t)


def xor(E,left,right):
    for a,b in zip(left,right):E.cx(a,b)


def peel_pair(m,E,math,W,inverse=False):
    a,u,*tail=W
    if inverse:
        if tail:xor(E,u,tail[0])
        scale(m,E,[u],a)
        if tail:math.inverse([u,*tail])
        scale(m,E,tail,u)
    else:
        scale(m,E,tail,u,True)
        if tail:math.inverse([u,*tail])
        scale(m,E,[u],a,True)
        if tail:xor(E,u,tail[0])


def extract(m,E,math,W,inverse=False):
    for k in (reversed(range(0,len(W),2)) if inverse else range(0,len(W),2)):
        peel_pair(m,E,math,W[k:],inverse)


def rational_layout(W):
    return W[::2],W[1::2]


def synthesize(m,E,W,inverse=False):
    if len(W)==2:return rational_layout(W)
    a,b=W[:2];P,Q=rational_layout(W[2:]);t=len(P)
    def fold(back=False):
        if back:
            scale(m,E,P,a,True);xor(E,P[0],b)
            for i in range(t-1):xor(E,P[i+1],Q[i])
            scale(m,E,Q,b,True)
            for p,q in zip(P,Q):xor(E,q,p)
        else:
            for p,q in zip(P,Q):xor(E,q,p)
            scale(m,E,Q,b)
            for i in range(t-1):xor(E,P[i+1],Q[i])
            xor(E,P[0],b);scale(m,E,P,a)
    if inverse:fold(True);synthesize(m,E,W[2:],True)
    else:synthesize(m,E,W[2:]);fold()
    return [a,*P],[b,*Q]


def model(F,s):
    """Independent field operations for the totalized circuit map."""
    W=list(s)
    for k in range(0,len(W),2):
        a,u=W[k:k+2];tail=W[k+2:];ui=F.inv(u or 1)
        normalized=[1]+[F.mul(x,ui) for x in tail]
        reciprocal=F.reciprocal(normalized,len(normalized))
        b=F.mul(u,F.inv(a or 1));out=reciprocal[1:]
        if out:out[0]^=b
        W[k+1]=b;W[k+2:]=out
    def build(W):
        if len(W)==2:return [W[0]],[W[1]]
        a,b=W[:2];P,Q=build(W[2:]);R=[p^q for p,q in zip(P,Q)]
        outQ=[F.mul(q,b or 1)^(R[i+1] if i+1<len(R) else 0) for i,q in enumerate(Q)]
        return [a]+[F.mul(r,a or 1) for r in R],[b^R[0]]+outQ
    P,Q=build(W)
    return trim([1]+Q),trim(P)


def emit(N,mode='both',count_only=False):
    m=load();m.inv_word=bind(m);b=m.Builder();W=[b.qubits(12) for _ in range(N)]
    class Counter(m.Em):
        extra_T=0
        def emit_mac(self,rt):self.extra_T+=sum(x[0]=='P' for x in rt)
        def emit_free(self,rt):pass
    cls=Counter if count_only else m.Em
    E=cls(b,*b.qubits(3),m.Alloc(b));math=BoundedArithmetic(m,E)
    def total():return b._cur_stats['toffoli']+getattr(E,'extra_T',0)
    stages=[]
    for name,fn in [('extract',lambda:extract(m,E,math,W)),('synthesize',lambda:synthesize(m,E,W)),
                    ('synthesize_inverse',lambda:synthesize(m,E,W,True)),('extract_inverse',lambda:extract(m,E,math,W,True))]:
        if mode=='synthesis' and name.startswith('extract'):continue
        before=total();fn();stages.append(dict(name=name,T=total()-before,peak=b.peak_qubits))
        if name=='synthesize':hook=len(b._items)
    report=dict(N=N,Q=b.peak_qubits,T=total(),stages=stages,
                public_Q=12*N,scratch_Q=b.peak_qubits-12*N,final_live=b.active_qubits,
                scope='exact arithmetic-emitter profile, including pivot tests, normalization, restoration; generic charts')
    return m,b,W,hook,report


def build(N,M,allocated_accumulators=False):
    m,b,W,hook,report=emit(N)
    tag=f'cfraction_n{N}_m{M}'+('_allocated_scan' if allocated_accumulators else '')
    stem=HERE/'candidates'/tag
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    tape=stem.with_suffix('.npy');np.save(tape,rows)
    kb=FastBuilder(km.CircuitBuilder())
    S=[kb.create_quantum_register(12,f'S{i}') for i in range(N)]
    masks=[kb.create_classical_register(12,f'E{i}') for i in range(M)]
    q=[x for w in S for x in w]+list(kb.alloc_qubits(b.next_qubit-12*N))
    fw,bk,tr=isomorphism(m)
    for word in S:apply_map(kb,word,map_ops(fw))
    kb.info_tag('CF_PHASE_BEGIN');kb.free(km.array(q[12*N:]))
    P,Q=rational_layout(W)
    P=[km.array([q[i] for i in word]) for word in P]
    Q=[km.array([q[i] for i in word]) for word in Q]
    gf=km.GF2Field(12)
    class Forney(RecursiveForney):gf=km.GF2Field(12)
    body=Forney()
    # Two independent locator functionals fit in B1/B2, and the numerator
    # fits in Omega0. This is a reversible linear change of coordinates.
    # It avoids three accumulator words in addition to the polynomial banks.
    if allocated_accumulators:
        # Preserve the first validated storage policy for exact reproduction.
        root=kb.alloc_qubits(12);derivative=kb.alloc_qubits(12);numerator=kb.alloc_qubits(12)
        def lin(word,target,c):linear_xor(kb,word,target,[tr(m.gmul(c,1<<i),bk) for i in range(12)])
        for j,mask in enumerate(masks):
            x=tr(gf.pow(2,(-j)%4095),fw)
            def evaluate():
                kb.x(root[0])
                for i,word in enumerate(Q,1):
                    lin(word,root,m.gpow(x,i))
                    if i%2:lin(word,derivative,m.gpow(x,i-1))
                for i,word in enumerate(P):lin(word,numerator,m.gpow(x,i))
            evaluate();body._forney_kernel(kb,root,derivative,numerator,mask);evaluate()
        kb.free(root);kb.free(derivative);kb.free(numerator)
    else:
        assert N>=4
        for word in S:apply_map(kb,word,map_ops(bk))
        derivative,root,numerator=Q[0],Q[1],P[0]
        def lin(word,target,c):linear_xor(kb,word,target,[gf.mul(c,1<<i) for i in range(12)])
        for j,mask in enumerate(masks):
            x=gf.pow(2,(-j)%4095);rootops=map_ops([gf.mul(gf.mul(x,x),1<<i) for i in range(12)])
            def evaluate(back=False):
                def d():
                    for i in range(3,len(Q)+1,2):lin(Q[i-1],derivative,gf.pow(x,i-1))
                def r():
                    lin(derivative,root,x)
                    for i in range(4,len(Q)+1,2):lin(Q[i-1],root,gf.pow(x,i))
                    kb.x(root[0])
                def p():
                    for i in range(1,len(P)):lin(P[i],numerator,gf.pow(x,i))
                if back:p();r();apply_map(kb,root,rootops,True);d()
                else:d();apply_map(kb,root,rootops);r();p()
            evaluate();body._forney_kernel(kb,root,derivative,numerator,mask);evaluate(True)
        for word in S:apply_map(kb,word,map_ops(fw))
    kb.info_tag('CF_PHASE_END')
    for word in S:apply_map(kb,word,map_ops(fw),True)
    wrapper=stem.with_suffix('.wrapper.kmx');kb.finish_circuit().to_file(str(wrapper),format='kmx')
    parts=[stem.with_suffix('.'+x+'.kmx') for x in ['pre','scan','post']];handles=[x.open('w') for x in parts];part=0
    for line in wrapper.read_text().splitlines(True):
        if 'CF_PHASE_BEGIN' in line:part=1
        elif 'CF_PHASE_END' in line:part=2
        else:handles[part].write(line)
    for h in handles:h.close()
    marks=HERE/'parameters/empty';marks.write_text('');out=stem.with_suffix('.kmx')
    subprocess.run([str(TRANSLATOR),str(tape),str(hook),str(12*M),*[str(p) for p in parts],str(out),str(b.next_out),str(marks),str(len(rows))],check=True,stdout=subprocess.DEVNULL)
    circ=km.Circuit.from_file(str(out),format='kmx');circ.to_file(str(stem.with_suffix('.kmb')),format='kmb')
    report.update(M=M,complete_Q=circ.num_qubits,complete_static_T=circ.operation_statistics().num_toffoli_like,
                  sha256=hashlib.sha256(stem.with_suffix('.kmb').read_bytes()).hexdigest())
    assert report['complete_static_T']==report['T']+80*M
    (HERE/'results'/f'{tag}_build.json').write_text(json.dumps(report,indent=2)+'\n')
    return circ,report


def test(N,M,shots,salt):
    circ,report=build(N,M);F=Field();seed=int.from_bytes(hashlib.sha256(salt.encode()).digest()[:8],'little')
    rng=np.random.default_rng(seed);S=[];truth=[];pred=[];charts=[];mask=rng.integers(0,4096,(shots,M))
    for i in range(shots):
        d=i%(N//2+1);pos=rng.choice(M,d,replace=False).tolist();a=[F.pow(2,j) for j in pos];e=rng.integers(1,4096,d).tolist()
        if d>=2 and i%11==0:
            z=0
            for aa,ee in zip(a[:-1],e[:-1]):z^=F.mul(aa,ee)
            v=F.mul(z,F.inv(a[-1]))
            if v:e[-1]=v
        s=F.syndrome(a,e,N);b,p=model(F,s);chart=exact_peel(F,s) is not None
        if chart:assert b==[1]+F.scale(F.locator(a),F.inv(F.locator(a)[0]))[1:] if d else b==[1]
        expected=0
        for j in range(M):
            beta=F.pow(2,(-j)%4095)
            if F.evaluate(b,beta)==0:
                denominator=F.evaluate(F.derivative(b),beta)
                value=F.mul(F.evaluate(p,beta),F.inv(denominator)) if denominator else 0
                expected^=(int(mask[i,j])&value).bit_count()%2
        S.append(s);pred.append(expected);truth.append(sum((int(mask[i,j])&v).bit_count() for j,v in zip(pos,e))%2);charts.append(chart)
    S=np.array(S,dtype=np.uint16)
    bits=lambda v:((np.asarray(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    public={q for reg in circ.register_data.values() for q in reg if isinstance(q,km.q)}
    scratch=km.array([km.q(i) for i in range(circ.num_qubits) if km.q(i) not in public])
    errors=dict(value=0,scratch=0,model_phase=0,mask=0,twice_phase=0,twice_value=0,twice_scratch=0)
    phases=[]
    for start in range(0,shots,256):
        stop=min(start+256,shots);sim=km.Simulator(batch_size=stop-start,seed=seed+start);sim.use_same_registers_as(circ)
        for i in range(N):sim.write_across_shots(circ.register_data[f'S{i}'],bits(S[start:stop,i]))
        for j in range(M):sim.write_across_shots(circ.register_data[f'E{j}'],bits(mask[start:stop,j]))
        sim.do(circ);phase=sim.read_phase_flipped_across_shots();phases.extend(phase.tolist())
        errors['model_phase']+=int(np.count_nonzero(phase^np.array(pred[start:stop],dtype=bool)))
        for i in range(N):errors['value']+=int(np.any(sim.read_across_shots(circ.register_data[f'S{i}'])^bits(S[start:stop,i]),axis=0).sum())
        for j in range(M):errors['mask']+=int(np.any(sim.read_across_shots(circ.register_data[f'E{j}'])^bits(mask[start:stop,j]),axis=0).sum())
        errors['scratch']+=int(np.any(sim.read_across_shots(scratch),axis=0).sum());sim.do(circ)
        errors['twice_phase']+=int(sim.read_phase_flipped_across_shots().sum())
        for i in range(N):errors['twice_value']+=int(np.any(sim.read_across_shots(circ.register_data[f'S{i}'])^bits(S[start:stop,i]),axis=0).sum())
        errors['twice_scratch']+=int(np.any(sim.read_across_shots(scratch),axis=0).sum())
    assert not any(errors.values()),errors
    truth_bad=np.array(phases)^np.array(truth,dtype=bool)
    assert not truth_bad[np.array(charts)].any()
    report.update(status='PASS_COMPONENT_CONTRACT',shots=shots,salt=salt,errors=errors,truth_phase_errors=int(truth_bad.sum()),
                  chart_exceptions=charts.count(False),truth_phase_errors_on_supported_charts=int(truth_bad[np.array(charts)].sum()),
                  scope='complete SMALL C-fraction oracle: actual host-basis syndromes, runtime masks, phase and restoration; forced unsupported gaps included; not benchmark acceptance')
    (HERE/f'results/cfraction_n{N}_m{M}_{salt}.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);p.add_argument('--m',type=int,default=31);p.add_argument('--shots',type=int,default=1024);p.add_argument('--salt',default='core-phase-cfraction-a');p.add_argument('--profile',action='store_true');a=p.parse_args()
    if a.profile:
        _,_,_,_,r=emit(a.n,count_only=True);(HERE/f'results/cfraction_n{a.n}_profile.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
    else:test(a.n,a.m,a.shots,a.salt)
