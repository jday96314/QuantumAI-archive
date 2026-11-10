"""Independent segmented-phase proof and measured conditional carry erasure.

Also quantify the incumbent's signed multirow carry width, exhaustively.
"""
import hashlib,json,math
from fractions import Fraction
from pathlib import Path
import numpy as np
import kickmix as km
from independent_arithmetic import controlled_add

def carry(b,a,t,out=None,phase=False):
 initial=b.alloc_qubits(1)[0]
 for i in range(len(a)):
  c=initial if i==0 else a[i-1]
  b.cx(a[i],t[i]);b.cx(a[i],c);b.ccx(c,t[i],a[i])
 if phase:b.z(a[-1])
 else:b.cx(a[-1],out)
 for i in reversed(range(len(a))):
  c=initial if i==0 else a[i-1]
  b.ccx(c,t[i],a[i]);b.cx(a[i],c);b.cx(a[i],t[i])
 b.free(km.array([initial]))

def kernel(n,split=None,measured=True):
 b=km.CircuitBuilder();a=b.create_quantum_register(n,'a');t=b.create_quantum_register(n,'t');one=b.alloc_qubits(1)[0];b.x(one)
 def apply(reg,offset=0):
  for i,q in enumerate(reg):b.z_pow(q,Fraction(2)**(i+offset-n-9))
 if split is None:
  out=b.alloc_qubits(n);b.cx(t,out);controlled_add(b,a,out,one);apply(out);controlled_add(b,a,out,one,True);b.cx(t,out);b.free(out)
 else:
  c=b.alloc_qubits(1)[0];low=b.alloc_qubits(split);b.cx(t[:split],low)
  carry(b,a[:split],t[:split],c)
  controlled_add(b,a[:split],low,one);apply(low);controlled_add(b,a[:split],low,one,True);b.cx(t[:split],low);b.free(low)
  high=b.alloc_qubits(n-split);b.cx(t[split:],high)
  controlled_add(b,a[split:],high,one,carry_in=c);apply(high,split);controlled_add(b,a[split:],high,one,True,carry_in=c)
  b.cx(t[split:],high);b.free(high)
  if measured:
   m=b.alloc_bits(1)[0];b.hmr(c,m);b.push_condition(m);carry(b,a[:split],t[:split],phase=True);b.pop_condition();b.bit_store(m,False);b.free(km.array([m]))
  else:carry(b,a[:split],t[:split],c)
  b.free(km.array([c]))
 b.x(one);b.free(km.array([one]));return b.finish_circuit()

def native():
 rows=[];rng=np.random.default_rng(731)
 for n,split in [(8,None),(8,4),(27,None),(27,14)]:
  c=kernel(n,split);N=8192;av=rng.integers(0,1<<n,N,dtype=np.uint64);tv=rng.integers(0,1<<n,N,dtype=np.uint64)
  s=km.Simulator(batch_size=N,counting_operations=True,seed=86753);s.use_same_registers_as(c);bits=lambda v:((v[None,:]>>np.arange(n,dtype=np.uint64)[:,None])&1).astype(bool)
  s.write_across_shots('a',bits(av));s.write_across_shots('t',bits(tv));s.do(c)
  assert np.array_equal(s.read_across_shots('a'),bits(av));assert np.array_equal(s.read_across_shots('t'),bits(tv))
  expected=((av+tv)&((1<<n)-1)).astype(object)
  theta=np.array([int(s.read_shot_phase(i)*(1<<63))&((1<<64)-1) for i in range(N)],dtype=object)
  wanted=np.array([int(v)<<(64-n-10) for v in expected],dtype=object)
  assert np.array_equal(theta,wanted),(n,split,theta[:10],wanted[:10])
  regs=set(c.register_data['a'])|set(c.register_data['t']);scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in regs]);assert not s.read_across_shots(scratch).any()
  row=dict(n=n,split=split,Q=c.num_qubits,static_Toffolis=c.operation_statistics().num_toffoli_like,mean_executed_Toffolis=s.operation_statistics().num_toffoli_like/N,shots=N,status='PASS');rows.append(row);print(json.dumps(row),flush=True)
 Path('results/segmented_native.json').write_text(json.dumps(rows,indent=2)+'\n')

def incumbent_carry():
 # Read-only reuse for incumbent probe, separate from independent seeds.
 import binary_centered_invsqrt as p
 table=np.load('candidates/margin_curve8_t5.9e-08.npy');cuts=[10,12,14,16,18,20];ranges={cut:[1<<60,-(1<<60)] for cut in cuts};residual=[1<<60,-(1<<60)];slope_range=[1<<60,-(1<<60)]
 def multiply(a,z,exponent):
  total=np.zeros_like(z);L=14;zz=z&((1<<L)-1)
  for g in range(7):
   bit=lambda k:0 if k<0 else (zz>>min(k,L-1))&1
   digit=bit(2*g-1)+bit(2*g)-2*bit(2*g+1)
   for mag in [1,2]:
    e=exponent+2*g+(mag==2);v=a<<e if e>=0 else (a>>-e)+((a>>(-e-1))&1)
    total+=np.where(abs(digit)==mag,np.sign(digit)*v,0)
  return total
 for start in range(1,3133453,131072):
  xs=np.arange(start,min(3133453,start+131072));k=21-(np.frexp(xs.astype(float))[1]-1);xx=xs<<k;addr=((xx>>13)&255)+((k&1)<<8);z=(xx&16383)-8192
  c0,c1,c2=table[addr].T;c1=c1-multiply(c2,z,-11);delta=-multiply(c1,z,-12)
  residual=[min(residual[0],int(delta.min())),max(residual[1],int(delta.max()))];slope_range=[min(slope_range[0],int(c1.min())),max(slope_range[1],int(c1.max()))]
  for cut in cuts:
   cc=((c0&((1<<cut)-1))+delta)>>cut;r=ranges[cut];ranges[cut]=[min(r[0],int(cc.min())),max(r[1],int(cc.max()))]
 rows=[]
 for cut,(lo,hi) in ranges.items():
  bits=next(n for n in range(1,32) if -(1<<(n-1))<=lo and hi<(1<<(n-1)))
  rows.append(dict(low_bits=cut,carry_range=[lo,hi],signed_carry_bits=bits,low_plus_carry_bits=cut+bits))
 result=dict(exhaustive_inputs=3133452,residual_range=residual,updated_slope_range=slope_range,rows=rows);Path('results/segmented_incumbent_carry.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':native();incumbent_carry()
