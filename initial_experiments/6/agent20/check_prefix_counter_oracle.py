"""Full small-oracle native check, with explicit approximate-chart errors."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
import kickmix as km
import prefix_counter_oracle as oracle
from algebra import Field
from check_algebra import peel
import prefix_counter_reference as reference
HERE=Path(__file__).resolve().parent

def test(N=14,M=127,pairs=3,through=3,suffix=2,shots=2048):
 started=time.time();c,r=oracle.build(N,M,pairs,through,suffix);F=Field();seed=2026102613;rng=np.random.default_rng(seed)
 masks=rng.integers(0,4096,(shots,M));syndromes=[];truth=[];predicted=[];supported=[]
 for i in range(shots):
  if i%7==0:
   s=rng.integers(0,4096,N).tolist()
   if i%14==0:s[::2]=[0]*(N//2)
   want=None;chart=False
  else:
   d=i%(N//2+1);pos=rng.choice(M,d,replace=False).tolist();a=[F.pow(2,j) for j in pos];e=rng.integers(1,4096,d).tolist()
   if d>=2 and i%11==0:
    z=0
    for aa,ee in zip(a[:-1],e[:-1]):z^=F.mul(aa,ee)
    v=F.mul(z,F.inv(a[-1]))
    if v:e[-1]=v
   s=F.syndrome(a,e,N);want=sum((int(masks[i,j])&v).bit_count() for j,v in zip(pos,e))%2;chart=peel(F,s) is not None
  f,g,regular=reference.model(F,s,pairs,through,suffix);syndromes.append(s);predicted.append(reference.phase(F,f,g,masks[i]));truth.append(want);supported.append(chart and regular)
 S=np.asarray(syndromes);bits=lambda v:((np.asarray(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool)
 public={q for reg in c.register_data.values() for q in reg if isinstance(q,km.q)};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);private_bits=km.array([km.b(i) for i in range(12*M,c.num_bits)]);errors=dict(value=0,mask=0,scratch=0,private_bits=0,model_phase=0,twice_value=0,twice_phase=0,twice_scratch=0,twice_private_bits=0);phases=[]
 for start in range(0,shots,256):
  stop=min(start+256,shots);sim=km.Simulator(batch_size=stop-start,seed=seed+start);sim.use_same_registers_as(c)
  for j in range(N):sim.write_across_shots(c.register_data[f'S{j}'],bits(S[start:stop,j]))
  for j in range(M):sim.write_across_shots(c.register_data[f'E{j}'],bits(masks[start:stop,j]))
  sim.do(c);phase=sim.read_phase_flipped_across_shots();phases.extend(phase.tolist());errors['model_phase']+=int((phase^np.asarray(predicted[start:stop],dtype=bool)).sum())
  for j in range(N):errors['value']+=int(np.any(sim.read_across_shots(c.register_data[f'S{j}'])^bits(S[start:stop,j]),axis=0).sum())
  for j in range(M):errors['mask']+=int(np.any(sim.read_across_shots(c.register_data[f'E{j}'])^bits(masks[start:stop,j]),axis=0).sum())
  errors['scratch']+=int(np.any(sim.read_across_shots(scratch),axis=0).sum());errors['private_bits']+=int(np.any(sim.read_across_shots(private_bits),axis=0).sum());sim.do(c);errors['twice_phase']+=int(sim.read_phase_flipped_across_shots().sum())
  for j in range(N):errors['twice_value']+=int(np.any(sim.read_across_shots(c.register_data[f'S{j}'])^bits(S[start:stop,j]),axis=0).sum())
  errors['twice_scratch']+=int(np.any(sim.read_across_shots(scratch),axis=0).sum())
  errors['twice_private_bits']+=int(np.any(sim.read_across_shots(private_bits),axis=0).sum())
 assert not any(errors.values()),errors
 truth_errors=[i for i in range(shots) if truth[i] is not None and bool(phases[i])!=bool(truth[i])];supported_errors=[i for i in truth_errors if supported[i]];assert not supported_errors,supported_errors[:5]
 r.update(status='PASS_SMALL_ORACLE_CONTRACT',shots=shots,seed=seed,errors=errors,truth_phase_errors=len(truth_errors),truth_phase_errors_on_supported_charts=len(supported_errors),supported_cases=sum(supported),arbitrary_syndrome_cases=sum(x is None for x in truth),seconds=time.time()-started,scope='full small oracle: independent totalized map, host-basis input, runtime masks, phase, scratch, inversion; approximate-chart truth errors counted; no benchmark acceptance')
 (HERE/'results'/f'prefix_counter_small_b{pairs}_t{through}_checks.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--shots',type=int,default=2048);p.add_argument('--pairs',type=int,default=3);p.add_argument('--through',type=int,default=3);a=p.parse_args();test(pairs=a.pairs,through=a.through,shots=a.shots)
