"""Squarefree cofactor padding using controlled affine transvections.

The underlying polynomial adapter is unchanged. Its controlled binary linear
maps are synthesized as rank-one transvections, avoiding a Toffoli for each
individual constant-field CNOT. Tests include arbitrary coefficient banks and
arbitrary flags, rather than only successful decoding paths.
"""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
import kickmix as km
from import_published import load
import deferred_euclid_padding as padding
from controlled_linear import apply as linear

ROOT=Path(__file__).resolve().parent


def apply(m,E,words,flags,inverse=False):
 old=padding.emit_ops;padding.emit_ops=lambda bb,ops,flag,inverse=False:linear(bb,ops,flag,inverse)
 try:return padding.direct_repad(m,E,words,flags,inverse)
 finally:padding.emit_ops=old


def build(N,which):
 m=load();b=m.Builder();words=[b.qubits(12) for _ in range(N)];flags=b.qubits(N//2);E=m.Em(b,*b.qubits(3),m.Alloc(b))
 for direction in ([False,True] if which=='both' else [True,False] if which=='inverse_both' else [which=='inverse']):apply(m,E,words,flags,direction)
 rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);stem=f'transvection_repad_n{N}_{which}';tape=ROOT/f'parameters/{stem}.npy';np.save(tape,rows);kb=km.CircuitBuilder()
 for i in range(N):kb.create_quantum_register(12,f'S{i}')
 kb.create_quantum_register(N//2,'flags');pre=ROOT/f'parameters/{stem}.pre.kmx';kb.finish_circuit().to_file(str(pre),format='kmx');empty=ROOT/f'parameters/{stem}.empty';empty.write_text('');out=ROOT/f'candidates/{stem}.kmx'
 subprocess.run([str(ROOT/'translate_prefix'),str(tape),str(len(rows)),'0',str(pre),str(empty),str(empty),str(out),str(b.next_out),str(empty),str(len(rows))],check=True)
 circuit=km.Circuit.from_file(str(out),format='kmx');circuit.to_file(str(out.with_suffix('.kmb')),format='kmb');return circuit


def reference(values,flags):
 m=load();N=len(values);n=N//2;_,(bp,b)=padding.frames(list(range(N)));out=values.copy()
 for k in range(1,n,2):
  constant=padding.QUADRATICS[(k-1)//2]
  for bank in [bp,b]:
   original=out[bank].copy()
   for i in range(max(0,n-k-2),len(bank)):
    delta=np.zeros(values.shape[1],dtype=np.int64)
    for shift,c in [(1,1),(2,constant)]:
     source=i+shift
     if source<len(bank):delta^=np.array([m.gmul(c,int(x)) for x in original[source]])
     elif source==len(bank):delta^=c
    out[bank[i]]^=delta*flags[k]
 return out


def main():
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=70);p.add_argument('--shots',type=int,default=4096);a=p.parse_args();start=time.time();rng=np.random.default_rng(715137);values=rng.integers(0,4096,(a.n,a.shots));flags=rng.integers(0,2,(a.n//2,a.shots)).astype(bool);want=reference(values,flags);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
 for which in ['forward','inverse','both','inverse_both']:
  c=build(a.n,which);sim=km.Simulator(batch_size=a.shots,seed=417);sim.use_same_registers_as(c);initial=want if which=='inverse' else values;expected=want if which=='forward' else values
  for i,row in enumerate(initial):sim.write_across_shots(c.register_data[f'S{i}'],bits(row))
  sim.write_across_shots(c.register_data['flags'],flags);sim.do(c);fail=[]
  for i,row in enumerate(expected):
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[f'S{i}'])!=bits(row)).any(0))
   if len(bad):fail.append(dict(word=i,shots=bad[:8].tolist()))
  assert np.array_equal(sim.read_across_shots(c.register_data['flags']),flags)
  public={q for r in c.register_data.values() for q in r};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());report=dict(which=which,status='PASS' if not fail and not len(junk) and not len(phase) else 'FAIL',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,failures=fail,scratch_errors=junk[:8].tolist(),phase_errors=phase[:8].tolist());reports.append(report);print(json.dumps(report),flush=True)
 report=dict(status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope='N70 squarefree cofactor padding, arbitrary coefficients and arbitrary flag states, independent polynomial model, inverse and fresh compositions, scratch and phase; no full decoder acceptance',N=a.n,shots=a.shots,circuits=reports,seconds=time.time()-start);Path(f'results/transvection_repad_n{a.n}_check.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
