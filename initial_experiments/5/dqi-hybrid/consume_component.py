"""Concrete reversible BM discrepancy overwrite, including zero diagonal.

Uses andigu/uncompute's Apache-2.0 field MAC, reversible swaps, and inverse
recovery. This is a COMPONENT, not a complete decoder. It illustrates the
precise extra information a zero locator coefficient requires retaining.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import kickmix as km
import numpy as np
from import_published import load

HERE=Path(__file__).resolve().parent
R4=HERE.parents[2]/'initial_experiments/4/codex'


def apply(m,E,A,S,park,inverse=False):
    """(x,0) -> (a*x + b, [a=0]*x), x=S[-1], a=A[-1].

    Other syndrome words and every locator word are preserved. The inverse
    maps every forward image back, including a=0, x!=0, and b!=0. The full
    gate map also remains reversible on arbitrary incoming park values, but
    the displayed formula and clean-workspace contract assume park starts 0.
    """
    assert len(A)==len(S) and len(A)>=2
    target=S[-1];pivot=A[-1]
    def add_other_terms():
        for coefficient,source in zip(A[:-1],S[:-1]):
            E.xor(m.prod_rt(target,coefficient,source))
    with m.eq_test(E,pivot,0) as zero:
        if not inverse:E.swaps(zero,list(zip(target,park)))
        E.cx(zero,pivot[0])
        with m.inv_word(E,pivot) as pi,m.opened(E,'consume_product',12) as product:
            if inverse:add_other_terms()
            m.mon_slot(E,target,pi if inverse else pivot,pivot if inverse else pi,product)
            if not inverse:add_other_terms()
        E.cx(zero,pivot[0])
        if inverse:E.swaps(zero,list(zip(target,park)))


def build(L,which):
    m=load();b=m.Builder()
    A=[b.qubits(12) for _ in range(L+1)]
    S=[b.qubits(12) for _ in range(L+1)]
    park=b.qubits(12)
    E=m.Em(b,*b.qubits(3),m.Alloc(b))
    hook=None
    if which in ['forward','both','phase']:
        with b.phase('consume/forward'):apply(m,E,A,S,park)
        hook=len(b._items)
    if which in ['inverse','both','phase']:
        with b.phase('consume/inverse'):apply(m,E,A,S,park,True)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9)
    stem=HERE/f'components/consume_l{L}_{which}'
    tape=stem.with_suffix('.npy');np.save(tape,rows)
    kb=km.CircuitBuilder()
    for key,words in [('A',A),('S',S),('P',[park])]:
        for i in range(len(words)):kb.create_quantum_register(12,f'{key}{i}')
    if which=='phase':kb.create_classical_register(12,'M')
    pre=stem.with_suffix('.pre.kmx');kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=HERE/'components/empty';empty.write_text('')
    phase=stem.with_suffix('.phase.kmx')
    phase.write_text(''.join(f'PUSH_CONDITION if b{i}\nZ q{q}\nPOP_CONDITION\n' for i,q in enumerate(S[-1])) if which=='phase' else '')
    marks=stem.with_suffix('.marks');marks.write_text('')
    out=stem.with_suffix('.kmx')
    subprocess.run([str(R4/'translate_prefix'),str(tape),str(hook if which=='phase' else len(rows)),
        str(12 if which=='phase' else 0),str(pre),str(phase),str(empty),str(out),str(b.next_out),str(marks),str(len(rows))],check=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(stem.with_suffix('.kmb')),format='kmb')
    stats=dict(which=which,L=L,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,
        public_quantum_bits=24*(L+1)+12,private_scratch=c.num_qubits-24*(L+1)-12,
        ops=len(c),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        circuit_sha256=hashlib.sha256(stem.with_suffix('.kmb').read_bytes()).hexdigest(),
        phase_summary=b.phase_summary())
    return c,stats


def field_table(m):
    log=np.array(m.LOG,dtype=np.int32);exp=np.array(m.EXP[:4095],dtype=np.uint16)
    table=exp[(log[:,None]+log[None,:])%4095]
    table[0,:]=0;table[:,0]=0
    return table


def reference(table,A,S):
    terms=table[A,S]
    target=np.bitwise_xor.reduce(terms,axis=1)
    park=np.where(A[:,-1]==0,S[:,-1],0).astype(np.uint16)
    out=S.copy();out[:,-1]=target
    return out,park


def check_native(L,shots,salt):
    """Independent value/phase checks of all four emitted native circuits.

    Exhaustive pivot values and exhaustive parked-input values precede random
    operands. All 102400 shots are used by default. This is not host DQI
    acceptance: these are component registers, not a complete decoder.
    """
    m=load();table=field_table(m)
    seed=int.from_bytes(hashlib.sha256(salt.encode()).digest()[:16],'little')
    rng=np.random.default_rng(seed)
    A=rng.integers(0,4096,(shots,L+1),dtype=np.uint16)
    S=rng.integers(0,4096,(shots,L+1),dtype=np.uint16)
    # Force every exceptional degree and every source value on a zero pivot.
    a=min(4096,shots);A[:a,-1]=np.arange(a,dtype=np.uint16)
    a=min(4096,max(0,shots-4096));A[4096:4096+a,-1]=0;S[4096:4096+a,-1]=np.arange(a,dtype=np.uint16)
    for degree in range(L+1):
        lo=8192+degree*64;hi=min(shots,lo+64)
        if lo<hi:A[lo:hi,degree+1:]=0
    A[-min(64,shots):]=0
    masks=rng.integers(0,4096,shots,dtype=np.uint16)
    out,park=reference(table,A,S)
    bits=lambda v:((v[:,None]>>np.arange(12))&1).astype(bool).T
    results=[]
    for which in ['forward','inverse','both','phase']:
        circuit,stats=build(L,which)
        errors=CounterErrors()
        native_t=0
        public=set(q for reg in circuit.register_data.values() for q in reg if isinstance(q,km.q))
        scratch=km.array([km.q(q) for q in range(circuit.num_qubits) if km.q(q) not in public])
        for start in range(0,shots,256):
            stop=min(shots,start+256);sl=slice(start,stop)
            sim=km.Simulator(batch_size=stop-start,seed=seed+start,counting_operations=True)
            sim.use_same_registers_as(circuit)
            initial=out if which=='inverse' else S
            for key,values in [('A',A),('S',initial)]:
                for i in range(L+1):sim.write_across_shots(circuit.register_data[f'{key}{i}'],bits(values[sl,i]))
            if which=='inverse':sim.write_across_shots(circuit.register_data['P0'],bits(park[sl]))
            if which=='phase':sim.write_across_shots(circuit.register_data['M'],bits(masks[sl]))
            sim.do(circuit)
            native_t+=sim.operation_statistics().num_toffoli_like
            bad=np.zeros(stop-start,dtype=bool)
            for key,values in [('A',A),('S',out if which=='forward' else S)]:
                for i in range(L+1):bad|=np.any(sim.read_across_shots(circuit.register_data[f'{key}{i}'])^bits(values[sl,i]),axis=0)
            expected_park=park[sl] if which=='forward' else np.zeros(stop-start,dtype=np.uint16)
            bad|=np.any(sim.read_across_shots(circuit.register_data['P0'])^bits(expected_park),axis=0)
            if which=='phase':bad|=np.any(sim.read_across_shots(circuit.register_data['M'])^bits(masks[sl]),axis=0)
            junk=np.any(sim.read_across_shots(scratch),axis=0)
            expected=np.array([int(v).bit_count()%2 for v in (out[sl,-1]&masks[sl])],dtype=bool) if which=='phase' else np.zeros(stop-start,dtype=bool)
            phase=sim.read_phase_flipped_across_shots()^expected
            errors.add(bad,junk,phase)
        stats.update(shots=shots,salt=salt,mean_native_T=native_t/shots,**errors.values)
        results.append(stats);print(json.dumps(stats),flush=True)
        assert not any(errors.values.values()),stats
    report=dict(status='PASS',scope='reversible discrepancy overwrite COMPONENT, not complete DQI decoder or task-cost win',
        L=L,shots_per_circuit=shots,salt=salt,circuits=results,
        zero_pivot_cases=int(np.count_nonzero(A[:,-1]==0)),all_pivot_values_covered=shots>=4096,
        every_parked_input_value_covered=shots>=8192)
    (HERE/f'results/consume_l{L}_{salt}.json').write_text(json.dumps(report,indent=2)+'\n')


class CounterErrors:
    def __init__(self):self.values=dict(value_errors=0,scratch_errors=0,phase_errors=0)
    def add(self,bad,junk,phase):
        self.values['value_errors']+=int(bad.sum())
        self.values['scratch_errors']+=int(junk.sum())
        self.values['phase_errors']+=int(phase.sum())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--l',type=int,default=34);p.add_argument('--shots',type=int,default=102400);p.add_argument('--salt',default='independent-a');a=p.parse_args()
    check_native(a.l,a.shots,a.salt)
