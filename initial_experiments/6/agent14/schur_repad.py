"""Count common x padding and replace paired factors with distinct quadratics."""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from fast_dqi import FastBuilder
from controlled_linear import apply as controlled
from quadratic_schur_model import F,node_maps


def count(m,E,B,index,flag,inverse=False):
 b=E.b.raw;width=len(index)
 for i in (reversed(range(len(B)+1)) if inverse else range(len(B)+1)):
  code=i+1
  def detection():
   with m.eq_test(E,index,0) as ze:
    if i==len(B):E.cx(ze,flag)
    else:
     with m.eq_test(E,B[i],0) as z:E.x(z);E.b.ccx(ze,z,flag);E.x(z)
  def store():
   for bit,q in enumerate(index):
    if code>>bit&1:E.cx(flag,q)
  def query():b.flip_if_equal(km.array([km.q(q) for q in index]),code,target=km.q(flag))
  if inverse:query();store();detection()
  else:detection();store();query()


def repad(m,E,A,B,index,flag,inverse=False):
 n=len(A);assert n%2 and len(B)==n;b=E.b.raw;_,_,_,_,tr=setup();_,kap,_,_=node_maps(n)
 for k in (reversed(range(1,n,2)) if inverse else range(1,n,2)):
  # index = common x multiplicity + 1.
  b.flip_if_less_than(km.array([km.q(q) for q in index]),n-k+1,target=km.q(flag));E.x(flag);ops=[];c=tr(kap[(k-1)//2],setup()[2])
  for bank,monic in [(A,False),(B,True)]:
   for i in range(max(0,n-k-2),n):
    for shift,value in [(1,1),(2,c)]:
     source=i+shift
     if source<n:
      for bit,q in enumerate(bank[source]):
       v=m.gmul(value,1<<bit)
       for j in range(12):
        if v>>j&1:ops.append(('cx',q,bank[i][j]))
     elif source==n and monic:
      for j in range(12):
       if value>>j&1:ops.append(('x',bank[i][j]))
  controlled(E.b,ops,flag,inverse)
  E.x(flag);b.flip_if_less_than(km.array([km.q(q) for q in index]),n-k+1,target=km.q(flag))


def build(n,kind='repad',inverse=False):
 m=setup()[0];b=FastBuilder(km.CircuitBuilder());A=[b.create_quantum_register(12,f'A{i}') for i in range(n)] if kind=='repad' else [];B=[b.create_quantum_register(12,f'B{i}') for i in range(n)];index=b.create_quantum_register((n+1).bit_length(),'index');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));flag=ab.qubits(1)[0]
 (repad(m,E,[[q.id for q in w] for w in A],[[q.id for q in w] for w in B],[q.id for q in index],flag,inverse) if kind=='repad' else count(m,E,[[q.id for q in w] for w in B],[q.id for q in index],flag,inverse));ab.release_recycle('repad_flag',[flag]);ab.release_recycle('repad_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test(n=7,shots=256):
 m,_,fw,bk,tr=setup();rng=np.random.default_rng(142087+n);_,kap,_,_=node_maps(n);aa=[];bb=[];pa=[];pb=[];cc=[]
 for case in range(shots):
  d=case%(n+1);c=n-d;A=rng.integers(0,4096,d).tolist();B=rng.integers(0,4096,d).tolist()+[1];P=[0,1] if c%2 else [1]
  for k in range(1,n,2):
   if k>=d:P=F.product(P,[kap[(k-1)//2],1,1])
  A0=([0]*c+A+[0]*n)[:n];B0=([0]*c+B+[0]*n)[:n];A1=(F.product(P,A)+[0]*n)[:n];B1=(F.product(P,B)+[0]*n)[:n];aa.append(A0);bb.append(B0);pa.append(A1);pb.append(B1);cc.append(c+1)
 bits=lambda v,w=12:((np.asarray(v)[None,:]>>np.arange(w)[:,None])&1).astype(bool);rows=[]
 for kind in ['count','repad']:
  cf=build(n,kind);ci=build(n,kind,True)
  for inverse,c in [(False,cf),(True,ci)]:
   sim=km.Simulator(batch_size=shots,seed=142087+inverse);sim.use_same_registers_as(c);srcA,srcB=(pa,pb) if kind=='repad' and inverse else (aa,bb);dstA,dstB=(pa,pb) if kind=='repad' and not inverse else (aa,bb)
   for name,values in [('A',srcA),('B',srcB)]:
    if kind=='count' and name=='A':continue
    for i in range(n):sim.write_across_shots(c.register_data[f'{name}{i}'],bits([tr(v[i],fw) for v in values]))
   if kind=='repad' or inverse:sim.write_across_shots(c.register_data['index'],bits(cc,len(c.register_data['index'])))
   sim.do(c);errors={}
   for name,values in [('A',dstA),('B',dstB)]:
    if kind=='count' and name=='A':continue
    for i in range(n):
     if not np.array_equal(sim.read_across_shots(c.register_data[f'{name}{i}']),bits([tr(v[i],fw) for v in values])):errors[f'{name}{i}']=True
   assert np.array_equal(sim.read_across_shots(c.register_data['index']),bits(np.zeros(shots,dtype=int) if kind=='count' and inverse else cc,len(c.register_data['index'])));public={q for reg in c.register_data.values() for q in reg};junk=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);r=dict(kind=kind,inverse=inverse,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,scratch=bool(sim.read_across_shots(junk).any()),phase=bool(sim.read_phase_flipped_across_shots().any()));rows.append(r);print(json.dumps(r),flush=True);assert not errors and not r['scratch'] and not r['phase']
 result=dict(status='PASS',n=n,shots=shots*4,rows=rows,scope='all common monomial padding lengths, arbitrary proper numerator and monic true denominator, quantum index, exact coefficients, source/index/scratch/phase');Path(f'results/schur_repad_n{n}_native.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=7);a=p.parse_args();test(a.n)
