"""Compact independent circuits, bind native samples, exact full-domain margins.
This is diagnostic development, not full-default native acceptance.
"""
import collections,hashlib,heapq,json,math,re
from pathlib import Path
import numpy as np
import kickmix as km
from native_sample import sample
from family_seeds import coordinates,eval_pipeline,MAX_X
from invsqrt.verify_submission_invsqrt import ProblemSpecInvSqrt,reference_values,rotation_t_count

def compact_local(c,path):
 c.to_file(str(path.with_suffix('.raw.kmx')),format='kmx');lines=path.with_suffix('.raw.kmx').read_text().splitlines()
 refs=lambda s:[] if s.startswith(('INFO_TAG','REGISTER','APPEND_TO_REGISTER')) else list(map(int,re.findall(r'(?<!\S)q(\d+)',s)))
 last={q:i for i,s in enumerate(lines) for q in refs(s)};mapping={q:q for q in range(22)};pool=[];nextwire=22;depth=0;out=[]
 for i,s in enumerate(lines):
  ids=refs(s)
  for q in ids:
   if q not in mapping:
    if pool:mapping[q]=heapq.heappop(pool)
    else:mapping[q]=nextwire;nextwire+=1
  out.append(re.sub(r'(?<!\S)q(\d+)',lambda m:'q'+str(mapping[int(m[1])]),s) if ids else s)
  reset=depth==0 and s.startswith(('R ','HMR '))
  for q in set(ids):
   if q>=22 and (reset or last[q]==i):heapq.heappush(pool,mapping.pop(q))
  if s.startswith('PUSH_CONDITION'):depth+=1
  if s.startswith('POP_CONDITION'):depth-=1
 path.with_suffix('.kmx').write_text('\n'.join(out)+'\n');d=km.Circuit.from_file(str(path.with_suffix('.kmx')));d.to_file(str(path),format='kmb');return d

def exact_model(frac,widths,kind,table):
 spec=ProblemSpecInvSqrt.for_n_p(10);margin=1<<64;worst=(0.,0);min_x=0
 for start in range(1,MAX_X+1,131072):
  xs=np.arange(start,min(MAX_X+1,start+131072),dtype=np.int64);k,u,addr,z=coordinates(xs,frac);raw=eval_pipeline(table,z,addr,widths,frac,kind)
  theta=np.left_shift(raw.astype(np.uint64),(64-widths[0]-10+k//2).astype(np.uint64));target,tol=reference_values(spec,xs.astype(np.uint64))
  error=np.abs((theta-target).view(np.int64)).view(np.uint64);margins=tol.astype(np.int64)-error.astype(np.int64);idx=margins.argmin()
  if int(margins[idx])<margin:margin=int(margins[idx]);min_x=int(xs[idx])
  ratios=error/tol;i=ratios.argmax()
  if ratios[i]>worst[0]:worst=(float(ratios[i]),int(xs[i]))
 return dict(exhaustive_inputs=MAX_X,exact_min_margin_units=margin,min_margin_x=min_x,worst_error_over_tolerance=worst[0],worst_x=worst[1],status='INTEGER_PASS' if margin>0 else 'INTEGER_FAIL')

def compare_sample_model(c,table,frac,widths,kind):
 rng=np.random.default_rng(20261025);xs=rng.integers(1,MAX_X+1,4096,dtype=np.int64);xs[:22]=1<<np.arange(22)
 x=c.register_data['x'];sim=km.Simulator(batch_size=len(xs),seed=281);sim.use_same_registers_as(c);sim.write_across_shots(x,((xs[None,:]>>np.arange(22)[:,None])&1).astype(bool));sim.do(c)
 k,u,a,z=coordinates(xs,frac);raw=eval_pipeline(table,z,a,widths,frac,kind)
 expected=np.left_shift(raw.astype(np.uint64),(64-widths[0]-10+k//2).astype(np.uint64));actual=np.array([int(sim.read_shot_phase(i)*(1<<63))&((1<<64)-1) for i in range(len(xs))],dtype=np.uint64)
 assert np.array_equal(expected,actual),(kind,xs[expected!=actual][:10]);return len(xs)

if __name__=='__main__':
 rows=[]
 for family,frac,widths,kind in [('serial',8,(27,18,7),'serial'),('linear',10,(28,18),'linear_serial'),('direct',7,(27,18,8),'direct_serial')]:
  original=km.Circuit.from_file('candidates/independent_'+family+'.kmb');c=compact_local(original,Path('candidates/independent_'+family+'_compact.kmb'));table=np.load('parameters/native_'+family+'.npy')
  match=compare_sample_model(c,table,frac,widths,kind);check=exact_model(frac,widths,kind,table);native=sample(c);margin=check['exact_min_margin_units'];rot=native['rotations'];eps=math.pi*margin/(1<<64)/rot
  # One direct T or more, derived from native statistics (sample driver
  # reports total equivalents with the sampled epsilon; read direct T here).
  stats=c.operation_statistics();tbase=stats.num_t_like
  eq=native['mean_Toffolis']+(tbase+rot*rotation_t_count(eps))/4
  row=dict(family=family,scope='full objective PROJECTED from native sample resources and exhaustive exact integer margin; not full-native verified',raw_Q=original.num_qubits,Q=c.num_qubits,model_match_shots=match,integer=check,native=native,direct_T=tbase,rotation_eps=eps,toffoli_equivalents=eq,projected_objective=c.num_qubits*eq)
  rows.append(row);print(json.dumps(row),flush=True);Path('results/family_final_checkpoints.json').write_text(json.dumps(rows,indent=2)+'\n')
