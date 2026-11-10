"""Native components for Wronskian transport and direct moment phase."""
from pathlib import Path
import argparse,hashlib,json,subprocess
import numpy as np
import kickmix as km
from import_published import load
from low_workspace_inverse import bind
from cfraction_decoder import scale,xor
from h_phase import binary_map
from algebra import Field
from check_algebra import wronskian_root

HERE=Path(__file__).resolve().parent


def moment(m,E,B,S,inverse=False):
    if inverse:
        for i in range(1,len(B)):E.xor(m.prod_rt(S[0],B[i-1],S[-i]))
        scale(m,E,[S[0]],B[-1],True)
    else:
        scale(m,E,[S[0]],B[-1])
        for i in range(1,len(B)):E.xor(m.prod_rt(S[0],B[i-1],S[-i]))


def transport(m,E,B,H,lam,q1,inverse=False):
    binary_map(m,E,lam,2048);binary_map(m,E,q1,2048)
    if not inverse:scale(m,E,H[:-1],lam)
    for w,target in zip(B,H):E.xor(m.prod_rt(target,q1,w))
    if inverse:scale(m,E,H[:-1],lam,True)
    binary_map(m,E,q1,2048,True);binary_map(m,E,lam,2048,True)


def build(k,kind):
    m=load();m.inv_word=bind(m);b=m.Builder()
    if kind=='moment':
        B=[b.qubits(12) for _ in range(k)];S=[b.qubits(12) for _ in range(k)];extras=[]
        registers=[(f'B{i}',w) for i,w in enumerate(B)]+[(f'S{i}',w) for i,w in enumerate(S)]
    else:
        B=[b.qubits(12) for _ in range(k+1)];S=[b.qubits(12) for _ in range(k+1)];extras=[b.qubits(12),b.qubits(12)]
        registers=[(f'B{i}',w) for i,w in enumerate(B)]+[(f'S{i}',w) for i,w in enumerate(S)]+[('lam',extras[0]),('q1',extras[1])]
    E=m.Em(b,*b.qubits(3),m.Alloc(b))
    fn=lambda back:moment(m,E,B,S,back) if kind=='moment' else transport(m,E,B,S,*extras,back)
    fn(False);hook=len(b._items);forward=b._cur_stats['toffoli'];fn(True)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);stem=HERE/f'candidates/{kind}_k{k}_phase'
    tape=stem.with_suffix('.npy');np.save(tape,rows);kb=km.CircuitBuilder()
    for name,word in registers:kb.create_quantum_register(len(word),name)
    kb.create_classical_register(12,'M');pre=stem.with_suffix('.pre.kmx');kb.finish_circuit().to_file(str(pre),format='kmx')
    phase=stem.with_suffix('.phase.kmx');target=S[0]
    phase.write_text(''.join(f'PUSH_CONDITION if b{i}\nZ q{q}\nPOP_CONDITION\n' for i,q in enumerate(target)))
    empty=HERE/'parameters/empty';empty.write_text('');out=stem.with_suffix('.kmx')
    subprocess.run([str(HERE/'../../4/codex/translate_prefix'),str(tape),str(hook),'12',str(pre),str(phase),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,stdout=subprocess.DEVNULL)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(stem.with_suffix('.kmb')),format='kmb')
    return c,dict(kind=kind,k=k,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,forward_T=forward,
                  public_Q=sum(len(w) for _,w in registers),sha256=hashlib.sha256(stem.with_suffix('.kmb').read_bytes()).hexdigest())


def test(k,kind,shots,salt):
    c,r=build(k,kind);F=Field(12,0x10eb);seed=int.from_bytes(hashlib.sha256(salt.encode()).digest()[:8],'little');rng=np.random.default_rng(seed)
    values={};want=[]
    if kind=='moment':
        B=rng.integers(0,4096,(shots,k));B[:,-1]=rng.integers(1,4096,shots);S=rng.integers(0,4096,(shots,k))
        for i in range(shots):
            v=F.mul(int(B[i,-1]),int(S[i,0]))
            for j in range(1,k):v^=F.mul(int(B[i,j-1]),int(S[i,-j]))
            want.append(v)
    else:
        B=rng.integers(0,4096,(shots,k+1));B[:,-1]=1;S=np.zeros((shots,k+1),dtype=int)
        lam=rng.integers(1,4096,shots);q1=rng.integers(0,4096,shots)
        for i in range(shots):
            A=rng.integers(0,4096,k).tolist();H=wronskian_root(F,A,B[i].tolist());S[i,:len(H)]=H
            q=[int(rng.integers(0,4096)),int(q1[i]),1]
            newB=F.add(F.product(q,B[i].tolist()),F.scale(A,int(lam[i])))
            newH=wronskian_root(F,B[i].tolist(),newB)
            want.append(newH[0] if newH else 0)
        values.update(lam=lam,q1=q1)
    for i in range(B.shape[1]):values[f'B{i}']=B[:,i]
    for i in range(S.shape[1]):values[f'S{i}']=S[:,i]
    masks=rng.integers(0,4096,shots);values['M']=masks
    bits=lambda a:((np.asarray(a)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    sim=km.Simulator(batch_size=shots,seed=seed);sim.use_same_registers_as(c)
    for key,arr in values.items():sim.write_across_shots(c.register_data[key],bits(arr))
    sim.do(c);expected=np.array([(int(a)&int(b)).bit_count()%2 for a,b in zip(want,masks)],dtype=bool)
    assert np.array_equal(sim.read_phase_flipped_across_shots(),expected)
    for key,arr in values.items():assert np.array_equal(sim.read_across_shots(c.register_data[key]),bits(arr)),key
    public={q for reg in c.register_data.values() for q in reg if isinstance(q,km.q)};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
    assert not sim.read_across_shots(scratch).any()
    r.update(status='PASS',shots=shots,salt=salt,scope='native arithmetic/phase sandwich, dynamic mask, all source and scratch restoration; component only')
    (HERE/f'results/{kind}_k{k}_{salt}.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--k',type=int,default=3);p.add_argument('--kind',choices=['moment','transport'],default='moment');p.add_argument('--shots',type=int,default=1024);p.add_argument('--salt',default='core-phase-recurrence-a');a=p.parse_args();test(a.k,a.kind,a.shots,a.salt)
