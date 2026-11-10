"""Concrete Bennett bridge from six Padé coordinates to transfer entries.

This simple schedule reconstructs the six moments in a temporary bank,
computes all matrix entries, and clears the temporary bank. It deliberately
preserves the Padé coordinates. Its meaning requires nonsingular second and
third minors, so only (1,1,1) and (2,1) are claimed. This is a resource/correctness
component, not a dispatcher. The width/cost failure of this bridge is not a
lower bound on bridges that derive entries in another representation.
"""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np,kickmix as km
from import_published import load
from low_workspace_inverse import bind
from lookahead_six_quantum import apply as six
from matrix_transfer import native_reference

ROOT=Path(__file__).resolve().parent

def apply(m,E,Q,M):
    # All operations into M are XORs, so the same schedule erases its output.
    def xor(x,y):
        for a,b in zip(x,y):E.cx(a,b)
    def mac(x,y,z):E.xor(m.prod_rt(z,x,y))
    with m.opened(E,'matrix_original_moments',72) as bank:
        F=[bank[i*12:(i+1)*12] for i in range(6)]
        for q,f in zip(Q,F):xor(q,f)
        six(m,E,*F,True)
        e,d,a,_,_,_=F;v1,v0,u0,u1,pi,v2=Q;A,U,C,V=M
        xor(d,A[0]);mac(e,u1,A[0]);xor(e,A[1])
        xor(u0,U[0]);xor(u1,U[1]);E.x(U[2][0])
        xor(a,C[0]);mac(d,v2,C[0]);mac(e,v1,C[0])
        xor(d,C[1]);mac(e,v2,C[1]);xor(e,C[2])
        xor(v0,V[0]);xor(v1,V[1]);xor(v2,V[2])
        six(m,E,*F)
        for q,f in zip(Q,F):xor(q,f)

def emit(which):
    m=load();m.inv_word=bind(m);b=m.Builder();Q=[b.qubits(12) for _ in range(6)]
    M=[[b.qubits(12) for _ in range(k)] for k in [2,3,3,3]]
    E=m.Em(b,*b.qubits(3),m.Alloc(b))
    for _ in range(2 if which=='both' else 1):apply(m,E,Q,M)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);stem='matrix_coefficients_'+which
    tape=ROOT/'parameters'/f'{stem}.npy';np.save(tape,rows)
    kb=km.CircuitBuilder()
    for prefix,n in [('Q',6),('M',11)]:
        for i in range(n):kb.create_quantum_register(12,f'{prefix}{i}')
    pre=ROOT/'parameters'/f'{stem}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=ROOT/'parameters'/'empty';empty.write_text('');out=ROOT/'candidates'/f'{stem}.kmx'
    subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c

def test(shots):
    rng=np.random.default_rng(2026102421);mul,inv,prod,transfer=native_reference();values=[];expected=[];patterns={}
    while len(values)<shots:
        F=rng.integers(0,4096,6).tolist()
        if len(values)%2:F[0]=0
        e,d,a,b,c,f=F;D=mul(d,d)^mul(e,a)
        if not D:continue
        u0=mul(mul(d,b)^mul(a,a),inv(D));u1=mul(mul(e,b)^mul(d,a),inv(D))
        p=c^mul(a,u0)^mul(b,u1)
        if not p:continue
        alpha=mul(f^mul(d,mul(u0,u0))^mul(b,mul(u1,u1)),inv(p));z=mul(p,inv(D))
        V=[mul(alpha,u0)^mul(z,d),u0^mul(alpha,u1)^mul(z,e),u1^alpha]
        Q=[V[1],V[0],u0,u1,inv(p),V[2]]
        M,pat=transfer(F);key=','.join(map(str,pat));patterns[key]=patterns.get(key,0)+1
        # Arbitrary destination values check XOR semantics and erasure.
        initial=rng.integers(0,4096,11).tolist();target=[x^y for x,y in zip(initial,sum(M,[]))]
        values.append(Q+initial);expected.append(Q+target)
    from published_dqi import isomorphism
    m=load();fw,_,tr=isomorphism(m)
    convert=lambda rows:np.array([[tr(int(v),fw) for v in row] for row in np.asarray(rows).T])
    values,expected=convert(values),convert(expected);names=[f'Q{i}' for i in range(6)]+[f'M{i}' for i in range(11)]
    bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
    for which in ['forward','inverse','both']:
        c=emit(which);sim=km.Simulator(batch_size=shots,seed=839121)
        sim.use_same_registers_as(c);source=expected if which=='inverse' else values;target=expected if which=='forward' else values
        for name,row in zip(names,source):sim.write_across_shots(c.register_data[name],bits(row))
        sim.do(c);errors={}
        for name,row in zip(names,target):
            bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(row)).any(0))
            if len(bad):errors[name]=bad[:8].tolist()
        public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        junk=int(sim.read_across_shots(scratch).any(0).sum());phase=int(sim.read_phase_flipped_across_shots().sum())
        reports.append(dict(direction=which,status='PASS' if not errors and not junk and not phase else 'FAIL',errors=errors,scratch_errors=junk,phase_errors=phase,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    report=dict(status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope=__doc__,shots=shots,patterns=patterns,circuits=reports)
    (ROOT/'results'/'matrix_coefficients.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--shots',type=int,default=4096);a=p.parse_args();test(a.shots)
