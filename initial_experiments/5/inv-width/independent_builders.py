"""Native circuits for independent seeds. Reuses only independent_arithmetic.

These are reviewable first circuits, not accepted submission candidates.
"""
import argparse,json,math,time
from fractions import Fraction
from pathlib import Path
import numpy as np
import kickmix as km
from independent_arithmetic import controlled_add,rounded_product
from family_seeds import fit,model_check
from native_sample import sample

def normalize(b,x,flags,reverse=False):
 for j in (range(5) if reverse else reversed(range(5))):
  step=1<<j;indices=range(22-step)
  if reverse:
   for i in indices:b.cswap(flags[j],x[i],x[i+step])
   b.flip_if_equal(x[-step:],0,target=flags[j])
  else:
   b.flip_if_equal(x[-step:],0,target=flags[j])
   for i in reversed(indices):b.cswap(flags[j],x[i],x[i+step])

def load(b,address,values,width):
 target=b.alloc_qubits(width);b.init_lookup(table=list(map(int,values)),address=address,target=target);return target

def erase(b,address,values,target):
 pad=b.alloc_qubits(max(0,(1<<(len(address)//2))-len(target)))
 b.del_lookup(table=list(map(int,values)),address=address,target=km.array([*target,*pad]));b.free(pad);b.free(target)

def phase(b,word,flags,denominator,signed=False,negative=False):
 """Phase word*2^s/2^denominator; signed extension then reversible shift."""
 head=b.alloc_qubits(10);full=km.array([*word,*head]);n=len(full)
 if signed:b.cx(word[-1],head)
 occupied=len(word);occs=[]
 for j in range(4):
  step=1<<j;occs.append(occupied)
  indices=range(n-step) if signed else range(min(occupied,n-step))
  for i in reversed(indices):b.cswap(flags[j+1],full[i],full[i+step])
  if signed:
   for i in range(step):b.ccx(flags[j+1],full[-1],full[i])
  occupied+=step
 for i,q in enumerate(full):
  angle=Fraction(2)**(i-denominator+1)
  if negative:angle=-angle
  if signed and i==n-1:angle=-angle
  b.z_pow(q,angle)
 for j in reversed(range(4)):
  step=1<<j
  if signed:
   for i in range(step):b.ccx(flags[j+1],full[-1],full[i])
  indices=range(n-step) if signed else range(min(occs[j],n-step))
  for i in indices:b.cswap(flags[j+1],full[i],full[i+step])
 if signed:b.cx(word[-1],head)
 b.free(head)

def table_load_all(b,address,table,sizes):
 packed=[sum(int(v)<<sum(sizes[:j]) for j,v in enumerate(row)) for row in table]
 allreg=load(b,address,packed,sum(sizes));regs=[];offset=0
 for size in sizes:regs.append(allreg[offset:offset+size]);offset+=size
 return allreg,regs,packed

def schoolbook(b,a,z,target,reverse=False,signed_coordinate=False):
 rows=list(enumerate(z))
 for j,ctrl in (reversed(rows) if reverse else rows):
  src=a[:min(len(a),len(target)-j)]
  if len(src):controlled_add(b,src,target[j:],ctrl,reverse^(signed_coordinate and j==len(z)-1))

def square_into(b,z,target,one,reverse=False):
 ops=[]
 for i in range(len(z)):
  ops.append((z[i],one,2*i,False))
  for j in range(i+1,len(z)):ops.append((z[i],z[j],i+j+1,j==len(z)-1))
 for src,ctrl,e,sub in (reversed(ops) if reverse else ops):
  controlled_add(b,km.array([src]),target[e:],ctrl,reverse^sub)

def build(frac,widths,kind,table,dirty=False):
 b=km.CircuitBuilder();x=b.create_quantum_register(22,'x');flags=b.alloc_qubits(5);normalize(b,x,flags)
 address=km.array([*x[21-frac:21],flags[0]]);L=21-frac;z=x[:L];sizes=[int(v).bit_length() for v in table.max(axis=0)]
 if kind=='direct_serial':
  # Each polynomial term gets its own lifetime, phases add exactly.
  c0=load(b,address,table[:,0],sizes[0]);phase(b,c0,flags,widths[0]+10);erase(b,address,table[:,0],c0)
  slope=load(b,address,table[:,1],sizes[1]);linear=b.alloc_qubits(21)
  b.x(z[-1]);rounded_product(b,slope,z,linear,widths[0]-widths[1]-21,borrow=x if dirty else ());b.x(z[-1])
  phase(b,linear,flags,widths[0]+10,signed=True,negative=True)
  b.x(z[-1]);rounded_product(b,slope,z,linear,widths[0]-widths[1]-21,True,borrow=x if dirty else ());b.x(z[-1]);b.free(linear);erase(b,address,table[:,1],slope)
  curve=load(b,address,table[:,2],sizes[2]);sq=b.alloc_qubits(2*L)
  b.x(z[-1]);square_into(b,z,sq,x[21]);b.x(z[-1])
  product=b.alloc_qubits(sizes[2]+2*L);schoolbook(b,curve,sq,product)
  shift=widths[2]+42-widths[0];assert shift>0
  high=product[shift:];roundbit=product[shift-1]
  controlled_add(b,km.array([x[21]]),high,roundbit)
  phase(b,high,flags,widths[0]+10)
  controlled_add(b,km.array([x[21]]),high,roundbit,True)
  schoolbook(b,curve,sq,product,True);b.free(product)
  b.x(z[-1]);square_into(b,z,sq,x[21],True);b.x(z[-1]);b.free(sq);erase(b,address,table[:,2],curve)
 else:
  allreg,regs,packed=table_load_all(b,address,table,sizes)
  # Full-width accumulator padding is part of this independent first schedule.
  extra=[b.alloc_qubits(max(0,widths[j]+int(j==0)-sizes[j])) for j in range(len(sizes))]
  regs=[km.array([*r,*p]) for r,p in zip(regs,extra)]
  b.x(z[-1])
  for j in reversed(range(len(widths)-1)):rounded_product(b,regs[j+1],z,regs[j],widths[j]-widths[j+1]-21,True,borrow=x if dirty else ())
  b.x(z[-1]);phase(b,regs[0],flags,widths[0]+10);b.x(z[-1])
  for j in range(len(widths)-1):rounded_product(b,regs[j+1],z,regs[j],widths[j]-widths[j+1]-21,borrow=x if dirty else ())
  b.x(z[-1])
  for p in extra:b.free(p)
  erase(b,address,packed,allreg)
 normalize(b,x,flags,True);b.free(flags);return b.finish_circuit()

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--family',choices=['serial','linear','direct'],required=True);p.add_argument('--dirty',action='store_true');a=p.parse_args()
 if a.family=='serial':frac,widths,kind=8,(27,18,7),'serial'
 elif a.family=='linear':frac,widths,kind=10,(28,18),'linear_serial'
 else:frac,widths,kind=7,(27,18,8),'direct_serial'
 table=fit(frac,widths,kind);check=model_check(table,frac,widths,kind);np.save('parameters/native_'+a.family+'.npy',table)
 c=build(frac,widths,kind,table,a.dirty);a.family+=('_dirty' if a.dirty else '');path='candidates/independent_'+a.family+'.kmb';c.to_file(path,format='kmb')
 result=dict(family=a.family,model=check,Q=c.num_qubits,static_Toffolis=c.operation_statistics().num_toffoli_like)
 try:result['native']=sample(c)
 except Exception as e:result['native']=dict(status='FAIL',error=str(e))
 Path('results/native_'+a.family+'.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
