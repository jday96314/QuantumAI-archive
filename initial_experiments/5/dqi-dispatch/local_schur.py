"""Six-word routed Schur head for 111,21,12,3, retaining regular inversions.

Three ordinary Schur inversions are shared. Degree gaps rotate their first
nonzero pivot forward; triangular corrections recover the moved zero coordinates
and prepare pivot-one dummy slots. No transfer matrix or copied original head
is needed inside this component. Classification flags are supplied and retained;
global pending-head construction/completion detection remain separate.
"""
import argparse,json,subprocess
from pathlib import Path
import numpy as np,kickmix as km
from import_published import load
from low_workspace_inverse import bind
from core_components import TrimmedArithmetic
from matrix_transfer import native_reference
ROOT=Path(__file__).resolve().parent

def apply(m,E,math,F,flags,inverse=False):
    first2,first3,second2=flags
    def permute(words,order,flag,back=False):
        current=list(range(len(words)));ops=[]
        for i,want in enumerate(order):
            j=current.index(want)
            if j!=i:ops.append((i,j));current[i],current[j]=current[j],current[i]
        for i,j in (reversed(ops) if back else ops):E.swaps(flag,list(zip(words[i],words[j])))
    def frobenius(flag,source,target,q):
        rt=[];m.lin_rows_terms(rt,math.rows(1,q),source,target)
        E.xor([('P',(flag,),0,t[1],t[2],t[3]) for t in rt])
    def product(flag,a,b,target,q=1):
        with m.opened(E,'head_mask',12) as masked:
            rt=[];m.lin_rows_terms(rt,math.rows(1,q),a,masked)
            rt=[('P',(flag,),0,t[1],t[2],t[3]) for t in rt]
            E.mac(rt);E.xor(m.prod_rt(target,masked,b));E.free(rt)
    def cube(flag,a,target):product(flag,a,a,target,2)
    def high5(flag):
        product(flag,F[1],F[3],F[5],2);product(flag,F[2],F[1],F[5],2);product(flag,F[1],F[1],F[5],4)
    def high4(flag):
        frobenius(flag,F[2],F[4],2);product(flag,F[1],F[2],F[4],2);frobenius(flag,F[1],F[4],4)
    def selected_scale(flag,word,factor,back=False):
        with m.opened(E,'head_selected_scale',12) as D:
            E.x(D[0]);rt=[('P',(flag,),0,(a,),0,(b,)) for a,b in zip(factor,D)]
            E.mac(rt);E.cx(flag,D[0])
            with m.inv_word(E,D) as di,m.opened(E,'head_scale_product',12) as temp:m.mon_slot(E,word,di if back else D,D if back else di,temp)
            E.cx(flag,D[0]);E.x(D[0]);E.free(rt)
    def schur(words,back):
        with m.inv_word(E,words[0]) as iv,m.opened(E,'head_schur_product',12) as temp:
            if not back:
                for word in words[1:]:m.mon_slot(E,word,iv,words[0],temp)
            math.inverse(words)
            if back:
                for word in words[1:]:m.mon_slot(E,word,words[0],iv,temp)
    def first_rotate(back):
        if back:
            permute(F,[2,3,4,5,0,1],first3,True);permute(F,[1,2,3,4,5,0],first2,True)
        else:
            permute(F,[1,2,3,4,5,0],first2);permute(F,[2,3,4,5,0,1],first3)
    def restore_first(back):
        with m.opened(E,'head_first_gap',1) as ff:
            E.cx(first2,ff[0]);E.cx(first3,ff[0])
            if back:
                selected_scale(first3,F[4],F[0],True);high4(first3)
                selected_scale(ff[0],F[5],F[0],True);high5(ff[0])
            else:
                high5(ff[0]);selected_scale(ff[0],F[5],F[0])
                high4(first3);selected_scale(first3,F[4],F[0])
            E.cx(first3,ff[0]);E.cx(first2,ff[0])
    def prepare_second(back):
        if back:
            cube(first3,F[3],F[5]);frobenius(first3,F[3],F[4],2);E.cx(first3,F[4][0]);E.cx(first3,F[2][0])
            permute(F[2:],[2,0,3,1],first3,True);permute(F[2:],[1,2,0,3],first2,True)
        else:
            permute(F[2:],[1,2,0,3],first2);permute(F[2:],[2,0,3,1],first3)
            E.cx(first3,F[2][0]);E.cx(first3,F[4][0]);frobenius(first3,F[3],F[4],2);cube(first3,F[3],F[5])
    def restore_second(back):
        def two():frobenius(first2,F[3],F[4],2)
        def three(flag):cube(flag,F[3],F[5])
        if back:
            selected_scale(second2,F[5],F[2],True);three(second2)
            selected_scale(first2,F[5],F[2],True);three(first2)
            selected_scale(first2,F[4],F[2],True);two()
        else:
            two();selected_scale(first2,F[4],F[2]);three(first2);selected_scale(first2,F[5],F[2]);three(second2);selected_scale(second2,F[5],F[2])
    def prepare_third(back):
        if back:
            E.cx(first2,F[4][0]);E.cx(second2,F[4][0]);permute(F[4:],[1,0],second2,True);permute(F[4:],[1,0],first2,True)
        else:
            permute(F[4:],[1,0],first2);permute(F[4:],[1,0],second2);E.cx(second2,F[4][0]);E.cx(first2,F[4][0])
    def slot_layout(back):
        # Quadratic: scale-only, then two MACs. Cubic: empty slot, then
        # beta/q2 preparation, then two MACs. Matches fixed_schedule.hpp.
        if back:
            E.cx(first3,F[1][0]);permute(F,[2,4,0,1,3,5],first3,True)
            E.cx(second2,F[3][0]);permute(F[3:5],[1,0],second2,True)
            E.cx(first2,F[1][0]);permute(F[1:3],[1,0],first2,True)
        else:
            permute(F[1:3],[1,0],first2);E.cx(first2,F[1][0])
            permute(F[3:5],[1,0],second2);E.cx(second2,F[3][0])
            permute(F,[2,4,0,1,3,5],first3);E.cx(first3,F[1][0])
    if inverse:
        slot_layout(True)
        permute(F[2:],[2,3,0,1],first2,True);schur(F[4:],True);prepare_third(True);restore_second(True)
        schur(F[2:],True);permute(F[2:],[1,2,3,0],second2,True);prepare_second(True);restore_first(True);schur(F,True);first_rotate(True)
    else:
        first_rotate(False);schur(F,False);restore_first(False);prepare_second(False);permute(F[2:],[1,2,3,0],second2)
        schur(F[2:],False);restore_second(False);prepare_third(False);schur(F[4:],False);permute(F[2:],[2,3,0,1],first2)
        slot_layout(False)

def emit(which,embedded=False):
    m=load();m.inv_word=bind(m);b=m.Builder();F=[b.qubits(12) for _ in range(6)];flags=b.qubits(3)
    if embedded:b.qubits(933-75)
    E=m.Em(b,*b.qubits(3),m.Alloc(b));math=TrimmedArithmetic(m,E);math.virtual_zero=True
    for back in ([False,True] if which=='both' else [True,False] if which=='inverse_both' else [which=='inverse']):apply(m,E,math,F,flags,back)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);stem='local_schur_'+which+('_embedded' if embedded else '')
    tape=ROOT/'parameters'/f'{stem}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
    for i in range(6):kb.create_quantum_register(12,f'F{i}')
    kb.create_quantum_register(3,'flags');pre=ROOT/'parameters'/f'{stem}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=ROOT/'parameters'/'empty';empty.write_text('');out=ROOT/'candidates'/f'{stem}.kmx'
    subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c

def test(shots):
    rng=np.random.default_rng(2026102519);mul,inv,prod,_=native_reference();sources=[];targets=[];flags=[];pats={}
    def add(a,b):
        a=list(a)+[0]*max(0,len(b)-len(a))
        for i,v in enumerate(b):a[i]^=v
        while a and not a[-1]:a.pop()
        return a
    def scale(a,x):return [mul(v,x) for v in a]
    def euclid(F):
        p=[0]*6+[1];r=F[::-1]
        while r and not r[-1]:r.pop()
        pattern=[];quotients=[];degree=0
        while degree<3:
            if not r:raise ValueError()
            g=len(p)-len(r)
            if degree+g>3:raise ValueError()
            beta=mul(r[-1],inv(p[-1]));a=p.copy();q=[0]*(g+1);ri=inv(r[-1])
            while len(a)>=len(r):
                j=len(a)-len(r);v=mul(a[-1],ri);q[j]=v
                for i,x in enumerate(r):a[i+j]^=mul(v,x)
                while a and not a[-1]:a.pop()
            q=scale(q,beta);quotients.append((beta,q));pattern.append(g);p,r=r,add(scale(p,beta),prod(q,r));degree+=g
        return pattern,quotients
    while len(sources)<shots:
        kind=len(sources)%4;F=rng.integers(0,4096,6).tolist()
        if kind==1:F[0]=0
        elif kind==2:F[0]=F[0] or 1;F[2]=mul(mul(F[1],F[1]),inv(F[0]))
        elif kind==3:F[0]=F[1]=0
        try:pat,qs=euclid(F)
        except ValueError:continue
        if pat!=[[1,1,1],[2,1],[1,2],[3]][kind]:continue
        if kind==0:out=sum(([be,q[0]] for be,q in qs),[])
        elif kind==1:out=[qs[0][0],0,qs[0][1][1],qs[0][1][0],qs[1][0],qs[1][1][0]]
        elif kind==2:out=[qs[0][0],qs[0][1][0],qs[1][0],0,qs[1][1][1],qs[1][1][0]]
        else:out=[1,0,qs[0][0],qs[0][1][2],qs[0][1][1],qs[0][1][0]]
        sources.append(F);targets.append(out);flags.append([kind==1,kind==3,kind==2]);key=str(pat);pats[key]=pats.get(key,0)+1
    from published_dqi import isomorphism
    m=load();fw,_,tr=isomorphism(m);convert=lambda x:np.array([[tr(int(v),fw) for v in row] for row in np.asarray(x).T]);sources,targets=convert(sources),convert(targets);fbits=np.array(flags,dtype=bool).T
    bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
    for which in ['forward','inverse','both','inverse_both']:
        c=emit(which);sim=km.Simulator(batch_size=shots,seed=913812);sim.use_same_registers_as(c);src=targets if which in ['inverse','inverse_both'] else sources;dst=targets if which in ['forward','inverse_both'] else sources
        for i,row in enumerate(src):sim.write_across_shots(c.register_data[f'F{i}'],bits(row))
        sim.write_across_shots(c.register_data['flags'],fbits);sim.do(c);errors={}
        for i,row in enumerate(dst):
            bad=np.flatnonzero((sim.read_across_shots(c.register_data[f'F{i}'])!=bits(row)).any(0))
            if len(bad):errors[i]=bad[:8].tolist()
        public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=int(sim.read_across_shots(scratch).any(0).sum());phase=int(sim.read_phase_flipped_across_shots().sum())
        reports.append(dict(direction=which,status='PASS' if not errors and not junk and not phase else 'FAIL',errors=errors,scratch_errors=junk,phase_errors=phase,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    report=dict(status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope=__doc__,shots=shots,patterns=pats,circuits=reports)
    (ROOT/'results'/'local_schur.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--shots',type=int,default=4096);a=p.parse_args();test(a.shots)
