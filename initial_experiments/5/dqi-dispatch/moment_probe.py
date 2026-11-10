"""Native resource/correctness probe for width-feasible conditional moment heads.

The fixed schedule has ordinary six-word heads, a carried one-slot prefix with
eight logical moments, and a carried two-slot cubic prefix with ten. Known
leading zeros reduce live moment words. This probes the separate-chart fallback;
it is not the desired shared single Möbius map or a complete dispatch block.
"""
import argparse,json,subprocess
from pathlib import Path
import numpy as np,kickmix as km
from import_published import load
from low_workspace_inverse import bind
from conditional_middle import ConditionalArithmetic,mobius
from matrix_transfer import native_reference
ROOT=Path(__file__).resolve().parent

def emit(k,n,which,embedded=False,count=False,leading_zero=True):
    m=load();m.inv_word=bind(m);b=m.Builder();A=[b.qubits(12) for _ in range(k+1)];D=[b.qubits(12) for _ in range(k)]
    p=b.qubits(12);T=([()] if leading_zero else [])+[b.qubits(12) for _ in range(n-int(leading_zero))];flag=b.qubits(1)[0]
    public=12*(2*k+n+2-int(leading_zero))+1
    if embedded:b.qubits(933-public)
    class Counter(m.Em):
        T=0
        def emit_mac(self,rt):self.T+=sum(t[0]=='P' for t in rt)
        def emit_free(self,rt):pass
    E=(Counter if count else m.Em)(b,*b.qubits(3),m.Alloc(b));math=ConditionalArithmetic(m,E,flag)
    for back in ([False,True] if which=='both' else [True,False] if which=='inverse_both' else [which=='inverse']):mobius(m,E,math,A,D,p,T,flag,back)
    if count:return dict(Q=b.peak_qubits,T=E.T+b._cur_stats['toffoli'])
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);stem=f'moment_probe_k{k}_n{n}_{which}'+('_embedded' if embedded else '')
    tape=ROOT/'parameters'/f'{stem}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
    for name,number in [('A',k+1),('D',k),('P',1),('T',n-1)]:
        for i in range(number):kb.create_quantum_register(12,f'{name}{i}')
    kb.create_quantum_register(1,'flag');pre=ROOT/'parameters'/f'{stem}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=ROOT/'parameters'/'empty';empty.write_text('');out=ROOT/'candidates'/f'{stem}.kmx'
    subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c

def test(k,n,shots):
    rng=np.random.default_rng(2026102529+k+n);mul,inv,_,_=native_reference();values=rng.integers(0,4096,(2*k+n+1,shots));flags=np.arange(shots)%2==1
    values[0]=values[k+1]=1;values[2*k+1]=rng.integers(1,4096,shots);want=values.copy()
    def sp(a,b,n):
        r=[0]*n
        for i,x in enumerate(a[:n]):
            for j,y in enumerate(b[:n-i]):r[i+j]^=mul(x,y)
        return r
    def sd(a,b,n):
        r=[0]*n;bi=inv(b[0]);a=a+[0]*n
        for i in range(n):
            z=a[i]
            for j in range(1,min(i+1,len(b))):z^=mul(b[j],r[i-j])
            r[i]=mul(z,bi)
        return r
    for s in range(shots):
        if not flags[s]:continue
        A=values[:k+1,s].tolist();D=values[k+1:2*k+1,s].tolist();p=int(values[2*k+1,s]);T=[0]+values[2*k+2:,s].tolist()
        out=sd(sp(A,T,n),[inv(p),0]+sp(D,T,n-2),n);assert out[0]==0;want[2*k+2:,s]=out[1:]
    from published_dqi import isomorphism
    m=load();fw,_,tr=isomorphism(m);convert=lambda x:np.array([[tr(int(v),fw) for v in row] for row in x]);values,want=convert(values),convert(want)
    names=[f'A{i}' for i in range(k+1)]+[f'D{i}' for i in range(k)]+['P0']+[f'T{i}' for i in range(n-1)];bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
    for which in ['forward','inverse','both','inverse_both']:
        c=emit(k,n,which);sim=km.Simulator(batch_size=shots,seed=827631);sim.use_same_registers_as(c);src=want if which in ['inverse','inverse_both'] else values;dst=want if which in ['forward','inverse_both'] else values
        for name,row in zip(names,src):sim.write_across_shots(c.register_data[name],bits(row))
        sim.write_across_shots(c.register_data['flag'],flags[None,:]);sim.do(c);errors={}
        for name,row in zip(names,dst):
            bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(row)).any(0))
            if len(bad):errors[name]=bad[:8].tolist()
        public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=int(sim.read_across_shots(scratch).any(0).sum());phase=int(sim.read_phase_flipped_across_shots().sum())
        reports.append(dict(direction=which,status='PASS' if not errors and not junk and not phase else 'FAIL',errors=errors,scratch_errors=junk,phase_errors=phase,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    report=dict(status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope=__doc__,k=k,n=n,shots=shots,circuits=reports,embedded=emit(k,n,'forward',True,True))
    (ROOT/'results'/f'moment_probe_k{k}_n{n}.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--k',type=int,default=19);p.add_argument('--n',type=int,default=8);p.add_argument('--shots',type=int,default=512);a=p.parse_args();test(a.k,a.n,a.shots)
