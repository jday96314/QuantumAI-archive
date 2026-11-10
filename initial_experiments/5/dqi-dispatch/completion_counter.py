"""Absorb a block's unary completion suffix into a six-bit history counter.

Promise: before any completion the incoming count is zero; after completion
every later local flag is one. With that promise the updated count regenerates
the flags, allowing their wires to be released. This does not assume that an
isolated zero pivot is completion. Gaps must first be correctly distinguished.
"""
import argparse,json,subprocess
from pathlib import Path
import numpy as np,kickmix as km
from import_published import load

ROOT=Path(__file__).resolve().parent

def increment(m,E,counter,flag,inverse=False):
    # Toggle high bits while their original lower bits remain intact. Each
    # equality ladder is erased before any of its operand bits are changed.
    for i in reversed(range(1,len(counter))):
        operands=[flag]+list(counter[:i]);value=1 if inverse else (1<<len(operands))-1
        with m.eq_test(E,operands,value) as eq:E.cx(eq,counter[i])
    E.cx(flag,counter[0])

def flags_from_count(m,E,counter,flags):
    for j,flag in enumerate(flags):
        threshold=len(flags)-j
        E.x(flag)
        for k in range(threshold):
            with m.eq_test(E,counter,k) as eq:E.cx(eq,flag)

def apply(m,E,counter,flags,inverse=False):
    if inverse:flags_from_count(m,E,counter,flags)
    for flag in (reversed(flags) if inverse else flags):increment(m,E,counter,flag,inverse)
    if not inverse:flags_from_count(m,E,counter,flags)

def emit(width,which):
    m=load();b=m.Builder();counter=b.qubits(6);flags=b.qubits(width);E=m.Em(b,*b.qubits(3),m.Alloc(b))
    for back in ([False,True] if which=='both' else [True,False] if which=='inverse_both' else [which=='inverse']):apply(m,E,counter,flags,back)
    rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);stem=f'completion_counter_b{width}_{which}'
    tape=ROOT/'parameters'/f'{stem}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
    kb.create_quantum_register(6,'count');kb.create_quantum_register(width,'flags')
    pre=ROOT/'parameters'/f'{stem}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx')
    empty=ROOT/'parameters'/'empty';empty.write_text('');out=ROOT/'candidates'/f'{stem}.kmx'
    subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
    c=km.Circuit.from_file(str(out),format='kmx');c.to_file(str(out.with_suffix('.kmb')),format='kmb');return c

def test(width,shots):
    rng=np.random.default_rng(2026102427+width)
    pairs=[(0,n) for n in range(width+1)]+[(n,width) for n in range(1,36-width)]
    selection=rng.integers(0,len(pairs),shots);old=np.array([pairs[i][0] for i in selection]);nflags=np.array([pairs[i][1] for i in selection])
    flags=np.arange(width)[:,None]>=width-nflags
    after=old+nflags;reports=[]
    bits=lambda v:((v[None,:]>>np.arange(6)[:,None])&1).astype(bool)
    for which in ['forward','inverse','both','inverse_both']:
        c=emit(width,which);sim=km.Simulator(batch_size=shots,seed=829171+width);sim.use_same_registers_as(c)
        from_after=which in ['inverse','inverse_both'];want_after=which in ['forward','inverse_both']
        sim.write_across_shots(c.register_data['count'],bits(after if from_after else old))
        sim.write_across_shots(c.register_data['flags'],np.zeros_like(flags) if from_after else flags)
        sim.do(c)
        errors=int(((sim.read_across_shots(c.register_data['count'])!=bits(after if want_after else old)).any(0)|(sim.read_across_shots(c.register_data['flags'])!=(np.zeros_like(flags) if want_after else flags)).any(0)).sum())
        public={q for reg in c.register_data.values() for q in reg};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=int(sim.read_across_shots(scratch).any(0).sum());phase=int(sim.read_phase_flipped_across_shots().sum())
        reports.append(dict(direction=which,status='PASS' if not errors and not junk and not phase else 'FAIL',value_errors=errors,scratch_errors=junk,phase_errors=phase,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
    report=dict(status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope=__doc__,width=width,shots=shots,promise_states=len(pairs),circuits=reports)
    (ROOT/'results'/f'completion_counter_b{width}.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--block',type=int,default=3);p.add_argument('--shots',type=int,default=4096);a=p.parse_args();test(a.block,a.shots)
