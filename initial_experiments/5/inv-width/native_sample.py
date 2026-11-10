"""Native diagnostic only; never substitutes for the full default host verifier."""
import argparse, hashlib, json, math, time
from pathlib import Path
import numpy as np
import kickmix as km
from invsqrt.verify_submission_invsqrt import (ProblemSpecInvSqrt,validate_registers,to_bits,from_bits,reference_values,half_turns_to_units,rotation_t_count)

MAX_X=3133452

def sample(circuit,count=8192,stage=False):
 spec=ProblemSpecInvSqrt.for_n_p(10)
 x,scratch=validate_registers(circuit,spec)
 rng=np.random.default_rng(20261024)
 boundaries=[1,2,3,4,MAX_X,958102]
 for k in range(22):
  for d in [-2,-1,0,1,2]:
   v=(1<<k)+d
   if 1<=v<=MAX_X:boundaries.append(v)
 # Uniform deterministic diagnostic sample plus table/normalization edges.
 for k in range(22):
  for i in range(257):
   v=((256+i)*(1<<k))//256
   if 1<=v<=MAX_X:
    boundaries.extend(vv for vv in [v-1,v,v+1] if 1<=vv<=MAX_X)
 count=max(count,len(boundaries));count=(count+255)//256*256
 xs=np.concatenate([np.array(boundaries,dtype=np.uint64),rng.integers(1,MAX_X+1,count-len(boundaries),dtype=np.uint64)])
 target,tol=reference_values(spec,xs)
 digest=hashlib.sha256(bytes.fromhex('00')).hexdigest() # overwritten with canonical host circuit digest
 from invsqrt.verify_submission_invsqrt import _circuit_sha256
 digest=_circuit_sha256(circuit)
 sim=km.Simulator(batch_size=256,debug_print=False,counting_operations=True,seed=int(digest,16));sim.use_same_registers_as(circuit)
 worst=0.;worst_x=0;margin=1<<64;start=time.monotonic();stage_counts={}
 segments=None
 if stage:
  import profile_incumbent
  path=Path('results/incumbent_tagged.kmx');lines=path.read_text().splitlines();r=profile_incumbent.analyze(path);spans=r['stage_spans']
  edges=[(name,a,b) for name,a,b in spans if b is not None]
  segments=[];previous=0
  for name,a,b in edges:
   if a>previous:segments.append(('setup',km.Circuit('\n'.join(lines[previous:a]))))
   segments.append((name,km.Circuit('\n'.join(lines[a:b]))));previous=b
  if previous<len(lines):segments.append(('setup',km.Circuit('\n'.join(lines[previous:]))))
 for offset in range(0,count,256):
  xx=xs[offset:offset+256];sim.clear_for_shot();sim.write_across_shots(x,to_bits(xx,22))
  if segments:
   for name,fragment in segments:
    before=sim.operation_statistics().num_toffoli_like;sim.do(fragment);delta=sim.operation_statistics().num_toffoli_like-before;stage_counts[name]=stage_counts.get(name,0)+delta
  else:sim.do(circuit)
  if not np.array_equal(from_bits(sim.read_across_shots(x)),xx):raise AssertionError('input not restored')
  dirty=sim.read_across_shots(scratch).any(axis=1)
  if dirty.any():raise AssertionError('dirty scratch '+str(scratch[int(np.argmax(dirty))]))
  theta=np.array([half_turns_to_units(sim.read_shot_phase(i)) for i in range(256)],dtype=np.uint64)
  used=np.abs((theta-target[offset:offset+256]).view(np.int64)).view(np.uint64);allowed=tol[offset:offset+256]
  ratios=used/allowed;idx=ratios.argmax()
  if ratios[idx]>worst:worst=float(ratios[idx]);worst_x=int(xx[idx])
  margin=min(margin,int((allowed.astype(object)-used.astype(object)).min()))
  if np.any(used>allowed):raise AssertionError(f'phase error x={xx[idx]} ratio={ratios[idx]}')
 stats=sim.operation_statistics();mean=stats.num_toffoli_like/count;tbase=stats.num_t_like/count;nrot=stats.num_z_pow(min_clifford_hierarchy_level=4)/count
 eps=math.pi*margin/(1<<64)/nrot if nrot else None;eq=mean+(tbase+nrot*rotation_t_count(eps))/4 if eps else mean+tbase/4
 return dict(status='SAMPLE_PASS',scope='diagnostic sampled-domain native; not accepted',Q=circuit.num_qubits,shots=count,mean_Toffolis=mean,static_Toffolis=circuit.operation_statistics().num_toffoli_like,rotations=nrot,sampled_min_margin_units=margin,sampled_worst_ratio=worst,worst_x=worst_x,sampled_objective=eq*circuit.num_qubits,toffoli_equivalent=eq,sha256=digest,seconds=time.monotonic()-start,stage_mean_Toffolis={k:v/count for k,v in stage_counts.items()})

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('circuit');p.add_argument('--count',type=int,default=8192);p.add_argument('--stages',action='store_true');p.add_argument('--output');a=p.parse_args()
 result=sample(km.Circuit.from_file(a.circuit),a.count,a.stages);print(json.dumps(result,indent=2))
 if a.output:Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
