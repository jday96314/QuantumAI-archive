"""Uniform reversible multiplication by a unimodular 2x2 polynomial matrix.

This component acts on equally sized low-coefficient bank interiors. It does
not graft consumed moments into the high boundary, move bank alignment, or
implement completion. It prices one concrete common arithmetic path, rather
than assuming quantum branch selection is free. V's monic cubic term is an
implicit constant. Matrix coefficient rows have lengths 2,3,3,3 (last V term
implicit), with nonzero constant determinant. All singular scalar pivots are
selected locally before division.
"""
from contextlib import contextmanager
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
from low_workspace_inverse import bind

ROOT=Path(__file__).resolve().parent

def apply(m,E,coeff,P,B,inverse=False):
    if not P:return
    A,U,C,V=coeff
    a,b,c,d=A[0],U[0],C[0],V[0]
    def mac(x,y,z): E.xor(m.prod_rt(z,x,y))
    def scale(w,x,xi,t): m.mon_slot(E,w,x,xi,t)
    with m.opened(E,'matrix_pivot',1) as ff:
        flag=ff[0]
        with m.eq_test(E,a,0) as z:E.cx(z,flag)
        # Only the constant matrix is factored; shifted terms keep their
        # original row order. Swap each output pair back before cross terms.
        E.swaps(flag,list(zip(a,c))+list(zip(b,d)))
        with m.inv_word(E,a) as ai,m.opened(E,'matrix_product',12) as t:
            scale(b,ai,a,t);mac(c,b,d);scale(c,ai,a,t)
            with m.inv_word(E,d) as di:
                def diagonal(i,back):
                    if back:
                        mac(c,P[i],B[i]);scale(B[i],di,d,t)
                        scale(P[i],ai,a,t);mac(b,B[i],P[i])
                    else:
                        mac(b,B[i],P[i]);scale(P[i],a,ai,t)
                        scale(B[i],d,di,t);mac(c,P[i],B[i])
                def cross(i):
                    for row,bank in [(A,P),(U,B)]:
                        for j,w in enumerate(row[1:],1):
                            if i>=j:mac(w,bank[i-j],P[i])
                    for row,bank in [(C,P),(V,B)]:
                        for j,w in enumerate(row[1:],1):
                            if i>=j:mac(w,bank[i-j],B[i])
                    if i>=3:
                        for q,tq in zip(B[i-3],B[i]):E.cx(q,tq)
                for i in (range(len(P)) if inverse else reversed(range(len(P)))):
                    if inverse:
                        cross(i);E.swaps(flag,list(zip(P[i],B[i])));diagonal(i,True)
                    else:
                        diagonal(i,False);E.swaps(flag,list(zip(P[i],B[i])));cross(i)
            scale(c,a,ai,t);mac(c,b,d);scale(b,a,ai,t)
        E.swaps(flag,list(zip(a,c))+list(zip(b,d)))
        with m.eq_test(E,a,0) as z:E.cx(z,flag)

def emit(n,which):
    m=load();m.inv_word=bind(m);b=m.Builder()
    rows=[[b.qubits(12) for _ in range(k)] for k in (2,3,3,3)]
    P=[b.qubits(12) for _ in range(n)];B=[b.qubits(12) for _ in range(n)]
    E=m.Em(b,*b.qubits(3),m.Alloc(b))
    directions=([False,True] if which=='both' else [True,False] if which=='inverse_both' else [which=='inverse'])
    for back in directions:apply(m,E,rows,P,B,back)
    tape=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    stem=f'matrix_transfer_n{n}_{which}';path=ROOT/'parameters'/f'{stem}.npy';np.save(path,tape)
    kb=km.CircuitBuilder()
    for group,count in [('M',11),('P',n),('B',n)]:
        for i in range(count):kb.create_quantum_register(12,f'{group}{i}')
    pre=ROOT/'parameters'/f'{stem}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=ROOT/'parameters'/'empty';empty.write_text('');out=ROOT/'candidates'/f'{stem}.kmx'
    subprocess.run([str(ROOT/'translate_prefix'),str(path),str(len(tape)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(tape))],check=True)
    circuit=km.Circuit.from_file(str(out),format='kmx');circuit.to_file(str(out.with_suffix('.kmb')),format='kmb')
    return circuit

def native_reference():
    ex=np.fromfile(ROOT/'../dqi-repair/parameters/antilog.bin',dtype=np.int32).tolist()
    lg=np.fromfile(ROOT/'../dqi-repair/parameters/log.bin',dtype=np.int32).tolist()
    def mul(a,b):return ex[(lg[a]+lg[b])%4095] if a and b else 0
    def inv(a):return ex[-lg[a]%4095] if a else 0
    def trim(a):
        while a and not a[-1]:a.pop()
        return a
    def add(a,b):
        a=list(a)+[0]*max(0,len(b)-len(a))
        for i,x in enumerate(b):a[i]^=x
        return trim(a)
    def prod(a,b):
        if not a or not b:return []
        r=[0]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):r[i+j]^=mul(x,y)
        return trim(r)
    def scaled(a,v):return trim([mul(x,v) for x in a])
    def divide(a,b):
        a=trim(list(a));q=[0]*max(0,len(a)-len(b)+1);bi=inv(b[-1])
        while len(a)>=len(b):
            k=len(a)-len(b);v=mul(a[-1],bi);q[k]=v
            for j,x in enumerate(b):a[k+j]^=mul(x,v)
            trim(a)
        return trim(q),a
    def transfer(F):
        p=[0]*6+[1];r=trim(F[::-1]);A,U,C,V=[1],[],[],[1];pattern=[]
        while sum(pattern)<3:
            if not r:raise ValueError('completion')
            q,rr=divide(p,r);g=len(q)-1
            if sum(pattern)+g>3:raise ValueError('straddles')
            beta=mul(r[-1],inv(p[-1]));q=scaled(q,beta)
            A,U,C,V=C,V,add(scaled(A,beta),prod(q,C)),add(scaled(U,beta),prod(q,V))
            p,r=r,add(scaled(p,beta),prod(q,r));pattern.append(g)
        det=add(prod(A,V),prod(U,C));assert len(det)==1 and det[0]
        assert V[-1]==1 and len(V)==4
        return [A+[0]*(2-len(A)),U+[0]*(3-len(U)),C+[0]*(3-len(C)),V[:3]],pattern
    return mul,inv,prod,transfer

def test(n,shots):
    begin=time.time();rng=np.random.default_rng(2026102417+n);mul,inv,prod,transfer=native_reference()
    values=[];wants=[];patterns={};kind_counts={}
    while len(values)<shots:
        F=rng.integers(0,4096,6).tolist();kind=len(values)%4
        if kind==1:F[0]=0;F[1]=F[1] or 1
        elif kind==2:
            F[0]=F[0] or 1;F[2]=mul(mul(F[1],F[1]),inv(F[0]))
        elif kind==3:F[0]=F[1]=0;F[2]=F[2] or 1
        try:M,pat=transfer(F)
        except ValueError:continue
        assert pat==[[1,1,1],[2,1],[1,2],[3]][kind],(kind,F,pat)
        key=','.join(map(str,pat));patterns[key]=patterns.get(key,0)+1
        pp=rng.integers(0,4096,n).tolist();bb=rng.integers(0,4096,n).tolist();A,U,C,V=M
        def total(a,b,c,d):
            x=prod(a,b)+[0]*n;y=prod(c,d)+[0]*n
            return [x[i]^y[i] for i in range(n)]
        np_=total(A,pp,U,bb);nb=total(C,pp,V+[1],bb)
        flat=sum(M,[]);values.append(flat+pp+bb);wants.append(flat+np_+nb)
    m=load();from published_dqi import isomorphism
    fw,_,tr=isomorphism(m)
    convert=lambda x:np.array([[tr(int(v),fw) for v in row] for row in np.asarray(x).T])
    values,wants=convert(values),convert(wants)
    bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    names=[f'M{i}' for i in range(11)]+[f'P{i}' for i in range(n)]+[f'B{i}' for i in range(n)]
    reports=[]
    for which in ['forward','inverse','both','inverse_both']:
        c=emit(n,which);sim=km.Simulator(batch_size=shots,seed=817128+n)
        sim.use_same_registers_as(c)
        source=wants if which in ['inverse','inverse_both'] else values
        target=wants if which in ['forward','inverse_both'] else values
        for name,row in zip(names,source):sim.write_across_shots(c.register_data[name],bits(row))
        sim.do(c);errors={}
        for name,row in zip(names,target):
            bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(row)).any(0))
            if len(bad):errors[name]=bad[:8].tolist()
        public={q for reg in c.register_data.values() for q in reg}
        scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        junk=int(sim.read_across_shots(scratch).any(0).sum());phase=int(sim.read_phase_flipped_across_shots().sum())
        reports.append(dict(direction=which,status='PASS' if not errors and not junk and not phase else 'FAIL',errors=errors,scratch_errors=junk,phase_errors=phase,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    report=dict(status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope=__doc__,n=n,shots=shots,patterns=patterns,circuits=reports,seconds=time.time()-begin)
    (ROOT/'results'/f'matrix_transfer_n{n}.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=9);p.add_argument('--shots',type=int,default=1024);a=p.parse_args();test(a.n,a.shots)
