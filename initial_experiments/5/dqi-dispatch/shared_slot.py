"""One shared field-product slot: scaling or an additive cofactor MAC.

For selector zero, X <- c X; for selector one, X <- X+cY, allowing c=0.
One selected-diagonal word is reused as the erasure MAC coefficient while
the diagonal inverse stays live, then restored before inverse cleanup.
This supplies two MACs in the final slot of a carried quadratic/cubic, without
allocating a second bank or matrix. Public c and Y are preserved.
"""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np,kickmix as km
from import_published import load
from low_workspace_inverse import bind
from monic_euclid_step import inverse_keep
from matrix_transfer import native_reference
ROOT=Path(__file__).resolve().parent

def conditional_inverse(m,E,word,mac_flag):
    with m.opened(E,'slot_conditional_inverse',12) as ci:
        inverse_keep(m,E,word,ci)
        E.x(mac_flag);E.swaps(mac_flag,list(zip(word,ci)));E.x(mac_flag)
        wins,rts=m.chain_open(E,word,ci);E.macs(rts[:-1]);E.free(rts[-1]);E.frees(rts[:-1]);m.chain_close(E,wins)

def apply(m,E,X,Y,c,pinv,mac_flag,inverse=False,norm=True):
    if not X and not norm:return
    if inverse and X:conditional_inverse(m,E,c,mac_flag)
    with m.opened(E,'slot_selected_diagonal',12) as diagonal:
        for a,b in zip(c,diagonal):E.cx(a,b);E.b.ccx(mac_flag,a,b)
        E.cx(mac_flag,diagonal[0])
        with m.inv_word(E,diagonal) as di,m.opened(E,'slot_product',12) as out:
            if norm:m.mon_slot(E,pinv,diagonal if inverse and not X else di,di if inverse and not X else diagonal,out)
            # Reuse diagonal as c on MAC inputs and zero on scaling inputs.
            for a,b in zip(c,diagonal):E.cx(a,b)
            E.cx(mac_flag,diagonal[0])
            for x,y in zip(X,Y):
                E.swaps(mac_flag,list(zip(x,y)));E.mac(m.prod_rt(out,c,x));E.swaps(mac_flag,list(zip(x,y)))
                E.mac([('P',(mac_flag,),0,(a,),0,(b,)) for a,b in zip(x,out)])
                E.free(m.prod_rt(x,di,out)+m.prod_rt(x,diagonal,y))
                E.swaps(None,list(zip(x,out)))
            E.cx(mac_flag,diagonal[0])
            for a,b in zip(c,diagonal):E.cx(a,b)
        E.cx(mac_flag,diagonal[0])
        for a,b in zip(c,diagonal):E.b.ccx(mac_flag,a,b);E.cx(a,b)
    if inverse and X:conditional_inverse(m,E,c,mac_flag)

def emit(n,which,embedded=False):
    m=load();m.inv_word=bind(m);b=m.Builder();X=[b.qubits(12) for _ in range(n)];Y=[b.qubits(12) for _ in range(n)];c=b.qubits(12);p=b.qubits(12);flag=b.qubits(1)[0]
    if embedded:b.qubits(933-(24*n+25))
    E=m.Em(b,*b.qubits(3),m.Alloc(b))
    for back in ([False,True] if which=='both' else [True,False] if which=='inverse_both' else [which=='inverse']):apply(m,E,X,Y,c,p,flag,back)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);stem=f'shared_slot_n{n}_{which}'+('_embedded' if embedded else '')
    tape=ROOT/'parameters'/f'{stem}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
    for name,k in [('X',n),('Y',n),('C',1),('P',1)]:
        for i in range(k):kb.create_quantum_register(12,f'{name}{i}')
    kb.create_quantum_register(1,'flag');pre=ROOT/'parameters'/f'{stem}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=ROOT/'parameters'/'empty';empty.write_text('');out=ROOT/'candidates'/f'{stem}.kmx'
    subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
    circuit=km.Circuit.from_file(str(out),format='kmx');circuit.to_file(str(out.with_suffix('.kmb')),format='kmb');return circuit

def test(n,shots):
    rng=np.random.default_rng(2026102513+n);mul,inv,_,_=native_reference();data=rng.integers(0,4096,(2*n+2,shots));flags=np.arange(shots)%2==1
    data[-2,~flags]=rng.integers(1,4096,(~flags).sum());data[-1]=rng.integers(1,4096,shots);data[-2,::7]=0;data[-2,(~flags)&(data[-2]==0)]=1
    want=data.copy()
    for s in range(shots):
        c=int(data[-2,s]);p=int(data[-1,s]);f=flags[s]
        for i in range(n):want[i,s]=int(data[i,s])^mul(c,int(data[n+i,s])) if f else mul(c,int(data[i,s]))
        want[-1,s]=p if f else mul(p,inv(c))
    from published_dqi import isomorphism
    m=load();fw,_,tr=isomorphism(m);convert=lambda a:np.array([[tr(int(v),fw) for v in row] for row in a]);data,want=convert(data),convert(want)
    names=[f'X{i}' for i in range(n)]+[f'Y{i}' for i in range(n)]+['C0','P0'];bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
    for which in ['forward','inverse','both','inverse_both']:
        c=emit(n,which);sim=km.Simulator(batch_size=shots,seed=917138);sim.use_same_registers_as(c)
        source=want if which in ['inverse','inverse_both'] else data;target=want if which in ['forward','inverse_both'] else data
        for name,row in zip(names,source):sim.write_across_shots(c.register_data[name],bits(row))
        sim.write_across_shots(c.register_data['flag'],flags[None,:]);sim.do(c);errors={}
        for name,row in zip(names,target):
            bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(row)).any(0))
            if len(bad):errors[name]=bad[:8].tolist()
        assert np.array_equal(sim.read_across_shots(c.register_data['flag']),flags[None,:])
        public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=int(sim.read_across_shots(scratch).any(0).sum());phase=int(sim.read_phase_flipped_across_shots().sum())
        reports.append(dict(direction=which,status='PASS' if not errors and not junk and not phase else 'FAIL',errors=errors,scratch_errors=junk,phase_errors=phase,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    report=dict(status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope=__doc__,n=n,shots=shots,circuits=reports)
    (ROOT/'results'/f'shared_slot_n{n}.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=9);p.add_argument('--shots',type=int,default=4096);a=p.parse_args();test(a.n,a.shots)
