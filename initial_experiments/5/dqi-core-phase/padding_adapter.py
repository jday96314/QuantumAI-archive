"""Explicit-leading cofactor repad and H unpad, quantum completion flags.

Adapts round-4 direct_repad and controlled_linear (both imported read-only).
Unlike its monic previous-cofactor bank, A/H have an explicit quantum lead.
"""
from pathlib import Path
import argparse,json,subprocess
import numpy as np
import kickmix as km
from import_published import load
from controlled_linear import apply as controlled
from quantum_coprime_padding import xor_constant,add_constant
from coprime_padding import QUADRATICS

HERE=Path(__file__).resolve().parent


def apply(m,E,A,B,flags,inverse=False):
    n=len(A)
    assert n%2 and len(flags)==n
    for k in (reversed(range(1,n,2)) if inverse else range(1,n,2)):
        ops=[];constant=QUADRATICS[(k-1)//2]
        for bank,monic in [(A,False)]+([(B,True)] if B is not None else []):
            for i in range(max(0,n-k-2),len(bank)):
                for shift,c in [(1,1),(2,constant)]:
                    src=i+shift
                    if src<len(bank):xor_constant(ops,c,bank[src],bank[i])
                    elif src==len(bank) and monic:add_constant(ops,c,bank[i])
        controlled(E.b,ops,flags[k],inverse)


def build(n,kind,which):
    m=load();b=m.Builder();A=[b.qubits(12) for _ in range(n)]
    B=[b.qubits(12) for _ in range(n)] if kind=='pair' else None
    flags=b.qubits(n);E=m.Em(b,*b.qubits(3),m.Alloc(b))
    for inverse in ([False,True] if which=='both' else [which=='inverse']):apply(m,E,A,B,flags,inverse)
    stem=HERE/f'candidates/padding_{kind}_n{n}_{which}';rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);tape=stem.with_suffix('.npy');np.save(tape,rows)
    kb=km.CircuitBuilder()
    for key,words in [('A',A)]+([('B',B)] if B is not None else []):
        for i in range(n):kb.create_quantum_register(12,f'{key}{i}')
    kb.create_quantum_register(n,'F');pre=stem.with_suffix('.pre.kmx');kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=HERE/'parameters/empty';empty.write_text('');out=stem.with_suffix('.kmx')
    subprocess.run([str(HERE/'../../4/codex/translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True,stdout=subprocess.DEVNULL)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(stem.with_suffix('.kmb')),format='kmb')
    return c,dict(n=n,kind=kind,which=which,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like)


def test(n,kind,shots):
    m=load();rng=np.random.default_rng(703415);keys=['A']+(['B'] if kind=='pair' else [])
    values={key:rng.integers(0,4096,(n,shots)) for key in keys};flags=rng.integers(0,2,(n,shots)).astype(bool)
    for i in range(min(n+1,shots)):flags[:,i]=np.arange(n)>=i
    expected={k:v.copy() for k,v in values.items()}
    for k in range(1,n,2):
        constant=QUADRATICS[(k-1)//2];table=np.array([m.gmul(constant,x) for x in range(4096)])
        for key,v in expected.items():
            old=v.copy()
            for i in range(max(0,n-k-2),n):
                delta=np.zeros(shots,dtype=np.int64)
                for shift,c in [(1,1),(2,constant)]:
                    src=i+shift
                    if src<n:delta^=old[src] if c==1 else table[old[src]]
                    elif src==n and key=='B':delta^=c
                v[i]^=delta*flags[k]
    bits=lambda a:((np.asarray(a)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
    reports=[]
    for which in ['forward','inverse','both']:
        c,r=build(n,kind,which);sim=km.Simulator(batch_size=shots,seed=79013);sim.use_same_registers_as(c)
        src=expected if which=='inverse' else values;want=expected if which=='forward' else values
        for key,v in src.items():
            for i,row in enumerate(v):sim.write_across_shots(c.register_data[f'{key}{i}'],bits(row))
        sim.write_across_shots(c.register_data['F'],flags);sim.do(c)
        for key,v in want.items():
            for i,row in enumerate(v):assert np.array_equal(sim.read_across_shots(c.register_data[f'{key}{i}']),bits(row)),(which,key,i)
        assert np.array_equal(sim.read_across_shots(c.register_data['F']),flags)
        public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public])
        assert not sim.read_across_shots(scratch).any();assert not sim.read_phase_flipped_across_shots().any()
        r.update(status='PASS',shots=shots,scope='quantum flags, explicit nonmonic A/H lead, arbitrary operands; forward, inverse, scratch and phase');reports.append(r);print(json.dumps(r),flush=True)
    (HERE/f'results/padding_{kind}_n{n}.json').write_text(json.dumps(reports,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=35);p.add_argument('--kind',choices=['pair','H'],default='pair');p.add_argument('--shots',type=int,default=512);a=p.parse_args();test(a.n,a.kind,a.shots)
