"""Shared in-place quotient/cofactor arithmetic, degrees 1, 2, and 3.

One target bank is scaled once and updated by quotient coefficients. Higher
coefficients use one reusable masked word. A cubic's x^3 term is an index-shift
of the same monomial path; the nonmonic partial cubic uses degree-two MACs.
No polynomial transfer matrix or Bennett bridge is formed. Previous quotient
gaps are represented by zero ghost words in the fixed target prefix.
"""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np,kickmix as km
from import_published import load
from low_workspace_inverse import bind
from matrix_transfer import native_reference
ROOT=Path(__file__).resolve().parent

def apply(m,E,A,B,Q,pinv,selectors,inverse=False,norm=True,max_degree=3):
    beta,q0,q1,q2=Q;g2,g3,p2,p3=selectors;d=len(B)
    assert len(A)==d-1 and d>=1
    def xor(a,b):
        for x,y in zip(a,b):E.cx(x,y)
    def mac(a,b,t):E.xor(m.prod_rt(t,a,b))
    def controlled_xor(flag,a,b):
        rt=[('P',(flag,),0,(x,),0,(y,)) for x,y in zip(a,b)];E.xor(rt)
    def original_lead():
        if d>=2:E.cx(p2,A[d-2][0])
        if d>=3:E.cx(p3,A[d-3][0])
    def remove_lead():
        if d>=2:controlled_xor(p2,A[d-2],beta)
        if d>=3:controlled_xor(p3,A[d-3],beta)
    def scale(back):
        with m.inv_word(E,beta) as bi,m.opened(E,'routed_scale_product',12) as temp:
            for a in A:m.mon_slot(E,a,bi if back else beta,beta if back else bi,temp)
            if norm:m.mon_slot(E,pinv,beta if back else bi,bi if back else beta,temp)
    def scalar():
        for i,a in enumerate(A):mac(q0,B[i],a)
        mac(q0,B[-1],beta)
    def monomial():
        targets=A+[beta,q0]
        for i,t in enumerate(targets):
            if i and i-1<d:xor(B[i-1],t)
        for flag,degree in [(g2,2),(g3,3)][:max_degree-1]:
            rt=[]
            for j,source in enumerate(B):
                destinations=[targets[i] for i in (j+1,j+degree) if i<len(targets)]
                for bit,x in enumerate(source):
                    out=tuple(w[bit] for w in destinations)
                    if out:rt.append(('P',(flag,),0,(x,),0,out))
            E.xor(rt)
    def higher(degree,back):
        if degree>max_degree:return
        word=q1 if degree==2 else q2
        with m.opened(E,'routed_masked_coefficient',12) as cache:
            def swap():
                with m.opened(E,'routed_selector',1) as ff:
                    E.cx(g3,ff[0])
                    if degree==2:E.cx(g2,ff[0])
                    E.swaps(ff[0],list(zip(cache,word)))
                    if degree==2:E.cx(g2,ff[0])
                    E.cx(g3,ff[0])
            def shear():
                if degree==2:
                    controlled_xor(g2,B[-1],cache)
                    if d>=2:controlled_xor(g3,B[-2],cache)
                else:controlled_xor(g3,B[-1],cache)
            def products():
                j=degree-1
                for i,a in enumerate(A):
                    if i>=j:mac(cache,B[i-j],a)
                if d-1-j>=0:mac(cache,B[d-1-j],beta)
                if d-j>=0:mac(cache,B[d-j],q0)
                if degree==3:mac(cache,B[-1],q1)
            swap()
            if back:shear()
            products()
            if not back:shear()
            swap()
    if inverse:
        monomial();higher(3,True);higher(2,True);scalar();remove_lead();scale(True);original_lead()
    else:
        original_lead();scale(False);remove_lead();scalar();higher(2,False);higher(3,False);monomial()
    return B,A+[beta,q0]+([q1] if max_degree>=2 else [])+([q2] if max_degree>=3 else [])

def emit(d,which,max_degree=3,embedded=False):
    m=load();m.inv_word=bind(m);b=m.Builder();A=[b.qubits(12) for _ in range(d-1)];B=[b.qubits(12) for _ in range(d)]
    Q=[b.qubits(12) for _ in range(4)];pi=b.qubits(12);flags=b.qubits(4)
    public=12*(2*d+4)+4
    if embedded:
        b.qubits(840+21+72-public)
    E=m.Em(b,*b.qubits(3),m.Alloc(b))
    for back in ([False,True] if which=='both' else [True,False] if which=='inverse_both' else [which=='inverse']):apply(m,E,A,B,Q,pi,flags,back,max_degree=max_degree)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);stem=f'routed_cofactor_d{d}_g{max_degree}_{which}'+('_embedded' if embedded else '')
    tape=ROOT/'parameters'/f'{stem}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
    for name,n in [('A',d-1),('B',d),('Q',4),('P',1)]:
        for i in range(n):kb.create_quantum_register(12,f'{name}{i}')
    kb.create_quantum_register(4,'flags')
    pre=ROOT/'parameters'/f'{stem}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=ROOT/'parameters'/'empty';empty.write_text('');out=ROOT/'candidates'/f'{stem}.kmx'
    subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c

def test(d,shots,max_degree=3):
    begin=time.time();rng=np.random.default_rng(2026102501+d+max_degree);mul,inv,prod,_=native_reference()
    inputs=[];outputs=[];flags=[];patterns={}
    for shot in range(shots):
        g=1+shot%max_degree;gp=1+(shot//max_degree)%3;h=d-gp
        bp=rng.integers(0,4096,h).tolist()+[1];bb=rng.integers(0,4096,d).tolist()+[1]
        aa=bp[:-1]+[0]*(d-1-h);q=rng.integers(0,4096,4).tolist();q[0]=int(rng.integers(1,4096));pi=int(rng.integers(1,4096))
        beta,q0,q1,q2=q;polynomial=[q0]+([q1] if g>=2 else [])+([q2] if g>=3 else [])+[1]
        result=prod(polynomial,bb)+[0]*(d+g+1)
        for i,v in enumerate(bp):result[i]^=mul(beta,v)
        result=result[:d+g+1];assert result[-1]==1
        out=aa+bb[:-1]+q+[pi]
        for i,v in enumerate(result[:d-1]):out[i]=v
        offset=2*d-1
        for j in range(g+1):out[offset+j]=result[d-1+j]
        out[-1]=mul(pi,inv(beta));inputs.append(aa+bb[:-1]+q+[pi]);outputs.append(out)
        flags.append([g==2,g==3,gp==2,gp==3]);key=f'{gp}->{g}';patterns[key]=patterns.get(key,0)+1
    from published_dqi import isomorphism
    m=load();fw,_,tr=isomorphism(m);convert=lambda x:np.array([[tr(int(v),fw) for v in row] for row in np.asarray(x).T])
    inputs,outputs=convert(inputs),convert(outputs);fbits=np.array(flags,dtype=bool).T
    bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    names=[f'A{i}' for i in range(d-1)]+[f'B{i}' for i in range(d)]+[f'Q{i}' for i in range(4)]+['P0'];reports=[]
    for which in ['forward','inverse','both','inverse_both']:
        c=emit(d,which,max_degree);sim=km.Simulator(batch_size=shots,seed=891721);sim.use_same_registers_as(c)
        src=outputs if which in ['inverse','inverse_both'] else inputs;dst=outputs if which in ['forward','inverse_both'] else inputs
        for name,row in zip(names,src):sim.write_across_shots(c.register_data[name],bits(row))
        sim.write_across_shots(c.register_data['flags'],fbits);sim.do(c);errors={}
        for name,row in zip(names,dst):
            bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(row)).any(0))
            if len(bad):errors[name]=bad[:8].tolist()
        assert np.array_equal(sim.read_across_shots(c.register_data['flags']),fbits)
        public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=int(sim.read_across_shots(scratch).any(0).sum());phase=int(sim.read_phase_flipped_across_shots().sum())
        reports.append(dict(direction=which,status='PASS' if not errors and not junk and not phase else 'FAIL',errors=errors,scratch_errors=junk,phase_errors=phase,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    report=dict(status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope=__doc__,degree=d,max_quotient=max_degree,shots=shots,patterns=patterns,circuits=reports,seconds=time.time()-begin)
    (ROOT/'results'/f'routed_cofactor_d{d}_g{max_degree}.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--degree',type=int,default=10);p.add_argument('--max-quotient',type=int,default=3);p.add_argument('--shots',type=int,default=4096);a=p.parse_args();test(a.degree,a.shots,a.max_quotient)
