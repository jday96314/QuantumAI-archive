"""Experimental cubic with compressed intermediate carries and borrowed phase headroom.

A consumed Horner source has its high carry integration undone before its
lookup-only high bits are measured. Replaying that small lookup and integration
restores the source for the inverse arithmetic. No native score is implied.
"""
import argparse,json
from fractions import Fraction
from pathlib import Path
import numpy as np
import kickmix as km
from binary_centered_invsqrt import cshift,del_flag_if_zero
from booth_adder import accumulate
from short_adder import add
from monomial_unlookup import erase
from optimize_workspace import compact
from phase_headroom import apply as phase


def lookup(b,address,target,values):
 values=list(map(int,values))
 def write(c,v):
  qs=km.array([q for i,q in enumerate(target) if v>>i&1])
  if len(qs):b.cx(c,qs)
 def visit(c,k,lo,hi):
  if all(v==values[lo] for v in values[lo:hi]):write(c,values[lo]);return
  bit=address[k];mid=(lo+hi)//2;tmp=b.alloc_qubits(1)[0];b.x(bit);b.init_and(c,bit,tmp);b.x(bit)
  visit(tmp,k-1,lo,mid);b.cx(c,tmp);visit(tmp,k-1,mid,hi);b.del_and(c,bit,tmp);b.free(km.array([tmp]))
 mid=len(values)//2;b.x(address[-1]);visit(address[-1],len(address)-2,0,mid);b.x(address[-1]);visit(address[-1],len(address)-2,mid,len(values))


def build(raw,limit=90,low=(19,11,4),product_limit=None,delay=2):
 f=6;lf=7;widths=(27,19,11,5);raw=raw.copy();index=np.arange(len(raw));parity=index>>lf;q=index&127
 correction=1<<17
 tab=raw.copy();tab[:,0]+=(q<<18)-(11<<23)-(parity<<25)-correction
 assert tab.min()>=0 and int(tab[:,0].max()).bit_length()<=23
 b=km.CircuitBuilder();x=b.create_quantum_register(22,'x');flags=b.alloc_qubits(5)
 def normalize(inv):
  for k in (range(5) if inv else reversed(range(5))):
   shift=1<<k
   if inv:cshift(b,x,shift,flags[k],22-shift,True);del_flag_if_zero(b,x[-shift:],flags[k])
   else:b.flip_if_equal(x[-shift:],0,target=flags[k]);cshift(b,x,shift,flags[k],22-shift)
 normalize(False);sizes=[int(c.max()).bit_length() for c in tab.T];regs=[b.alloc_qubits(s) for s in sizes];cutoff=sizes[0]-delay;address=km.array([*x[14:21],flags[0]])
 loadregs=[regs[0][:cutoff],*regs[1:]];loadvalues=tab.copy();loadvalues[:,0]&=(1<<cutoff)-1
 packed=[sum(int(v)<<sum(map(len,loadregs[:j])) for j,v in enumerate(row)) for row in loadvalues];lookup(b,address,km.array([q for r in loadregs for q in r]),packed)
 counters=[None]*3;head=km.array([]);pending={};oldhigh={};delayed_meas=None
 def spare():return km.array([*head,*regs[0][cutoff:]]) if delayed_meas is not None or not c0_loaded[0] else head
 c0_loaded=[False]
 def integrate(j,inv=False):
  high=regs[j][low[j]:];b.cx(x[14],high);add(b,counters[j],high,subtract=inv,carry_workspace=spare(),qubit_limit=limit);b.cx(x[14],high)
 def arith(j,inv=False):
  target=km.array([*regs[j][:low[j]],*counters[j]]);b.x(x[14]);b.x(x[21]);accumulate(b,regs[j+1],x[:15],target,widths[j]-widths[j+1]-21,subtract=not inv,rounded=True,borrow=spare(),auxiliary=km.array([x[21]]),qubit_limit=limit,product_limit=product_limit);b.x(x[21]);b.x(x[14])
 def drop(j,start):
  high=regs[j][start:];out=b.alloc_bits(len(high));b.hmr(high,out);pending[j]=out;oldhigh[j]=high;b.free(high);regs[j]=regs[j][:start]
 def restore(j,start):
  high=b.alloc_qubits(sizes[j]-start);regs[j]=km.array([*regs[j],*high]);lookup(b,address,high,tab[:,j]>>start)
  for cb,qb in zip(pending[j],high):b.push_condition(cb);b.z(qb);b.pop_condition()
  b.bit_store(pending[j],False);b.free(pending.pop(j))
 for j in [2,1,0]:
  counters[j]=b.alloc_qubits(1);arith(j)
  if j==0:lookup(b,address,regs[0][cutoff:],tab[:,0]>>cutoff);c0_loaded[0]=True
  integrate(j)
  if j==2:drop(3,0);head=b.alloc_qubits(3)
  if j==1:integrate(2,True);drop(2,low[2])
 y=km.array([*regs[0],*head,x[21],*x[:10]])
 assert len(y)==37
 def affine(inv=False):
  def constant(back):
   one=b.alloc_qubits(1);b.x(one)
   for bit in (reversed([17,23,24]) if back else [17,23,24]):add(b,one,y[bit:27],subtract=back,qubit_limit=limit)
   b.x(one);b.free(one)
  def parity(back):add(b,km.array([flags[0]]),y[25:27],subtract=back,qubit_limit=limit)
  def coarse(back):add(b,x[14:21],y[18:27],subtract=not back,qubit_limit=limit)
  for step in (reversed([constant,parity,coarse]) if inv else [constant,parity,coarse]):step(inv)
 affine();occ=[27,28,30,34]
 for k in range(4):cshift(b,y,1<<k,flags[k+1],occ[k])
 phase(b,y,27,flags[1:])
 for k in reversed(range(4)):cshift(b,y,1<<k,flags[k+1],occ[k],True)
 affine(True)
 for j in [0,1,2]:
  if j==1:restore(2,low[2]);integrate(2)
  if j==2:b.free(head);head=km.array([]);restore(3,0)
  integrate(j,True)
  if j==0:delayed_meas=b.alloc_bits(delay);b.hmr(regs[0][cutoff:],delayed_meas);c0_loaded[0]=False
  arith(j,True);b.free(counters[j]);counters[j]=None
 target=km.array([q for r in [regs[0][:cutoff],*regs[1:]] for q in r]);packed=[v|(int(row[0])>>cutoff)<<len(target) for v,row in zip(packed,tab)]
 erase(b,address,target,packed,extra_outcomes=delayed_meas);normalize(True)
 for reg in [*regs,flags]:b.free(reg)
 return b.finish_circuit()


def main():
 p=argparse.ArgumentParser();p.add_argument('--table',default='parameters/cubic6f7_w27_19_11_5_safe.txt');p.add_argument('--limit',type=int,default=94);p.add_argument('--output',default='candidates/cubic_phase_corrected_q94.kmb');p.add_argument('--product',type=int);a=p.parse_args();raw=np.loadtxt(a.table,skiprows=1,dtype=np.int64);path=Path(a.output);c=compact(build(raw,a.limit,product_limit=a.product),path);c.to_file(str(path),format='kmb');print(json.dumps(dict(status='BUILT',scope='requires quantum and full native verification',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,path=str(path))),flush=True)
if __name__=='__main__':main()
