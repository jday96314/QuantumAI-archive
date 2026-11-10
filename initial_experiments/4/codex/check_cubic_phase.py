"""Independent integer-model phase check for the experimental cubic circuit."""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
import binary_centered_invsqrt as poly
from invsqrt.verify_submission_invsqrt import half_turns_to_units


def main():
 p=argparse.ArgumentParser();p.add_argument('--candidate',default='candidates/cubic_phase_corrected_q94.kmb');p.add_argument('--shots',type=int,default=32768);p.add_argument('--table',default='parameters/cubic6f7_w27_19_11_5_safe.txt');a=p.parse_args();c=km.Circuit.from_file(a.candidate,format='kmb');rng=np.random.default_rng(95241);xs=np.r_[np.arange(1,min(a.shots,4096)+1),rng.integers(1,3133453,max(0,a.shots-4096))];n=len(xs)
 s=km.Simulator(batch_size=n,seed=712);s.use_same_registers_as(c);reg=c.register_data['x'];bits=((xs[None,:]>>np.arange(22)[:,None])&1).astype(bool);s.write_across_shots(reg,bits);s.do(c)
 source=np.flatnonzero(np.any(s.read_across_shots(reg)!=bits,0));scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in set(reg)]);junk=np.flatnonzero(s.read_across_shots(scratch).any(0));raw=np.loadtxt(a.table,skiprows=1,dtype=np.int64);poly.table_for=lambda *unused:raw;v=poly.model(xs,6,(27,19,11,5),rounded=True,lookup_frac=7,booth=True);want=v.astype(np.uint64)<<np.uint64(27);got=np.array([half_turns_to_units(s.read_shot_phase(i)) for i in range(n)],dtype=np.uint64);bad=np.flatnonzero(got!=want)
 out=dict(status='PASS' if not len(source) and not len(junk) and not len(bad) else 'FAIL',scope='cubic scratch/input/phase against independent integer model, no native acceptance',candidate=a.candidate,shots=n,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,source_errors=len(source),scratch_errors=len(junk),phase_errors=len(bad),examples=[dict(x=int(xs[j]),expected=int(want[j]),actual=int(got[j])) for j in bad[:8]])
 Path('results/'+Path(a.candidate).stem+'_quantum_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
