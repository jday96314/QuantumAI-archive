"""Prepare a short monotone table by binary search and coherent comparisons."""
import numpy as np
import kickmix as km
from monomial_unlookup import erase

def qrom(b,address,target,table):
 work=b.alloc_qubits(max(0,len(address)-1))
 def visit(c,level,lo,hi):
  if len(set(table[lo:hi]))==1:b.gf2_iadd(int(table[lo]),target=target,control=c);return
  bit=address[level];mid=(lo+hi)//2;t=work[level]
  b.x(bit);b.init_and(c,bit,t);b.x(bit);visit(t,level-1,lo,mid);b.cx(c,t);visit(t,level-1,mid,hi);b.del_and(c,bit,t)
 mid=len(table)//2;b.x(address[-1]);visit(address[-1],len(address)-2,0,mid);b.x(address[-1]);visit(address[-1],len(address)-2,mid,len(table));b.free(work)

def less_than(b,a,threshold,target,one):
 """target ^= a<threshold; one starts and ends at one. No scratch words."""
 n=len(a);assert len(threshold)==n
 b.x(threshold);b.x(target)
 for i in range(n-1):
  c=one if i==0 else threshold[i-1]
  b.cx(threshold[i],a[i]);b.cx(threshold[i],c);b.ccx(c,a[i],threshold[i])
 c=one if n==1 else threshold[n-2];last=threshold[-1]
 b.cx(last,a[-1]);b.cx(last,c);b.ccx(c,a[-1],target);b.cx(last,target);b.cx(last,c);b.cx(last,a[-1])
 for i in reversed(range(n-1)):
  c=one if i==0 else threshold[i-1]
  b.ccx(c,a[i],threshold[i]);b.cx(threshold[i],c);b.cx(threshold[i],a[i])
 b.x(threshold)

def prepare(b,address,target,table,workspace,one):
 """The table has two nonincreasing halves, selected by address[-1]."""
 n=len(address)-1;d=len(target);size=1<<n
 values=np.asarray(table,dtype=np.int64).reshape(2,size)
 assert np.all(np.diff(values,axis=1)<=0)
 assert len(workspace)>=n+1
 threshold=workspace[:n+1]
 for k in reversed(range(d)):
  prefix=target[k+1:];query=km.array([*prefix,address[-1]])
  lookup=[]
  for parity in range(2):
   for pre in range(1<<len(prefix)):
    trial=(pre<<(k+1))+(1<<k);lookup.append(int(np.sum(values[parity]>=trial)))
  qrom(b,query,threshold,lookup)
  less_than(b,address[:n],threshold[:n],target[k],one);b.cx(threshold[n],target[k])
  erase(b,query,threshold,lookup)

def test():
 import json
 rng=np.random.default_rng(8118);b=km.CircuitBuilder();a=b.create_quantum_register(9,'a');t=b.create_quantum_register(4,'t');w=b.alloc_qubits(9);one=b.alloc_qubits(1)[0];b.x(one)
 values=np.sort(rng.integers(0,16,(2,256)),axis=1)[:,::-1].flatten().tolist();prepare(b,a,t,values,w,one);b.x(one);b.free(km.array([one]));b.free(w);c=b.finish_circuit();N=1024
 aa=np.arange(N)%512;bits=lambda a,n:((a[None,:]>>np.arange(n)[:,None])&1).astype(bool)
 s=km.Simulator(batch_size=N,seed=516);s.use_same_registers_as(c);s.write_across_shots(a,bits(aa,9));s.do(c)
 assert np.array_equal(s.read_across_shots(a),bits(aa,9));assert np.array_equal(s.read_across_shots(t),bits(np.array(values)[aa],4));assert not s.read_phase_flipped_across_shots().any()
 public=set(a)|set(t);scratch=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);assert not s.read_across_shots(scratch).any()
 r=dict(status='PASS',shots=N,addresses=512,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like);open('results/monotone_lookup_check.json','w').write(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
if __name__=='__main__':test()
