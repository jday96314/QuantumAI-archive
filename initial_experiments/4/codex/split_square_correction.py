"""Add a rounded small square while keeping its high-bit source unchanged.

R is split at resbits-squarebits. Accumulate the monomial-rounded correction
of floor(R/2^split)^2 into only the low band and a signed three-bit counter.
Counter integration can then expose the complete corrected residual, and its
inverse makes the original source available again for counter cleanup.
"""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from short_adder import add
from goldschmidt_row_model import rounded,monomial_square


def apply(b,R,counter,split,shift,inverse=False,qubit_limit=None):
 high=R[split:];target=km.array([*R[:split],*counter]);steps=[]
 for i in range(len(high)):
  c=int(rounded(np.int64(3),shift-2*i))
  if c:steps.append((i,None,c))
  for j in range(i+1,len(high)):
   c=int(rounded(np.int64(3),shift-i-j-1))
   if c:steps.append((i,j,-c if j==len(high)-1 else c))
 for i,j,c in (reversed(steps) if inverse else steps):
  scratch=()
  if j is None:control=high[i]
  else:
   scratch=b.alloc_qubits(1);control=scratch[0];b.init_and(high[i],high[j],control)
  bits=[k for k in range(abs(c).bit_length()) if abs(c)>>k&1]
  if len(bits)==2 and bits[1]==bits[0]+1 and (qubit_limit is None or b.num_allocated_qubits+max(0,len(target)-bits[0]-2)<=qubit_limit):
   # The plain adder restores each source immediately after a carry step,
   # so repeated read-only source wires are valid. Avoid its in-place
   # prefix fallback here: that would temporarily modify aliased sources.
   lane=target[bits[0]:]
   lhs=km.array([control]*min(2,len(lane)))
   if len(lhs):add(b,lhs,lane,subtract=(c<0)^inverse)
  else:
   for k in (reversed(bits) if inverse else bits):
    if k<len(target):add(b,km.array([control]),target[k:],subtract=(c<0)^inverse,qubit_limit=qubit_limit)
  if len(scratch):b.del_and(high[i],high[j],control);b.free(scratch)


def integrate(b,R,counter,split,inverse=False,qubit_limit=None):
 high=R[split:]
 def unsigned(undo):add(b,counter,high,subtract=undo,qubit_limit=qubit_limit)
 def sign(undo):
  if len(high)>len(counter):add(b,km.array([counter[-1]]),high[len(counter):],subtract=not undo,qubit_limit=qubit_limit)
 for f in ([sign,unsigned] if inverse else [unsigned,sign]):f(inverse)


def build(bits,split,shift,which,ambient=0,limit=None):
 b=km.CircuitBuilder();R=b.create_quantum_register(bits,'R');c=b.create_quantum_register(3,'counter')
 if ambient:b.create_quantum_register(ambient,'ambient')
 for inv in ([False,True] if which=='both' else [True,False] if which=='inverse_both' else [which=='inverse']):
  if inv:integrate(b,R,c,split,True,limit);apply(b,R,c,split,shift,True,limit)
  else:apply(b,R,c,split,shift,False,limit);integrate(b,R,c,split,False,limit)
 return b.finish_circuit()


def main():
 p=argparse.ArgumentParser();p.add_argument('--bits',type=int,default=21);p.add_argument('--split',type=int,default=10);p.add_argument('--shift',type=int,default=10);p.add_argument('--shots',type=int,default=8192);p.add_argument('--ambient',type=int,default=0);p.add_argument('--limit',type=int);a=p.parse_args();rng=np.random.default_rng(741518);rr=rng.integers(0,1<<a.bits,a.shots);cc=rng.integers(0,8,a.shots);high=rr>>a.split;th=np.where(high>>(a.bits-a.split-1),high-(1<<(a.bits-a.split)),high);correction=monomial_square(th,a.bits-a.split,a.shift)
 low=((cc<<a.split)+(rr&((1<<a.split)-1))+correction)% (1<<(a.split+3));cout=low>>a.split;signed=np.where(cout&4,cout-8,cout);rout=(((high+signed)<<(a.split))+(low&((1<<a.split)-1)))% (1<<a.bits)
 bits=lambda v,n:((v[None,:]>>np.arange(n)[:,None])&1).astype(bool);reports=[]
 for which in ['forward','inverse','both','inverse_both']:
  c=build(a.bits,a.split,a.shift,which,a.ambient,a.limit);sim=km.Simulator(batch_size=a.shots,seed=741);sim.use_same_registers_as(c);initial=[rout,cout] if which=='inverse' else [rr,cc];want=[rout,cout] if which=='forward' else [rr,cc]
  for name,v in zip(['R','counter'],initial):sim.write_across_shots(c.register_data[name],bits(v,len(c.register_data[name])))
  sim.do(c);fail={}
  for name,v in zip(['R','counter'],want):
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[name])!=bits(v,len(c.register_data[name]))).any(0))
   if len(bad):fail[name]=bad[:8].tolist()
  public={q for r in c.register_data.values() for q in r};scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());r=dict(which=which,status='PASS' if not fail and not len(junk) and not len(phase) else 'FAIL',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,failures=fail,scratch_errors=junk[:8].tolist(),phase_errors=phase[:8].tolist());reports.append(r);print(json.dumps(r),flush=True)
 report=dict(status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',scope='split monomial square and signed carry integration, arbitrary incoming counter, inverse and fresh compositions, scratch and phase; no full inverse-square-root acceptance',parameters=vars(a),correction_range=[int(correction.min()),int(correction.max())],circuits=reports);tag='split_square_correction'+('_ambient' if a.ambient else '');Path('results/'+tag+'_check.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
