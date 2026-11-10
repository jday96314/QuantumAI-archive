"""Transpose the no-repair index with a completed-prefix repair index."""
import kickmix as km
from extension_inverse import Adapter,setup
from fast_dqi import FastBuilder


def build(N,inverse=False):
 m=setup()[0];n=N//2;b=FastBuilder(km.CircuitBuilder());S=[b.create_quantum_register(12,f'S{i}') for i in range(N)];index=b.create_quantum_register(max(1,(n-1).bit_length()),'index');ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab))
 for i in (reversed(range(n-1)) if inverse else range(n-1)):
  code=i+1;p=(code&-code).bit_length()-1
  for j,q in enumerate(index):
   if j!=p and code>>j&1:b.cx(index[p],q)
  other=[q.id for j,q in enumerate(index) if j!=p]
  with m.eq_test(E,[q.id for q in S[2*i+1]],0) as zc:
   if i:
    with m.eq_test(E,[q.id for q in S[2*i-1]],0) as zp,m.opened(E,'completion_boundary',1) as flag:
     E.x(zp);rt=[('P',(zc,),0,(zp,),0,(flag[0],))];E.mac(rt);E.x(zp)
     if len(other)>1:
      with m.eq_test(E,other,0) as zi:E.b.ccx(flag[0],zi,index[p].id)
     elif other:E.x(other[0]);E.b.ccx(flag[0],other[0],index[p].id);E.x(other[0])
     else:E.cx(flag[0],index[p].id)
     E.x(zp);E.free(rt);E.x(zp)
   else:
    if len(other)>1:
     with m.eq_test(E,other,0) as zi:E.b.ccx(zc,zi,index[p].id)
    elif other:E.x(other[0]);E.b.ccx(zc,other[0],index[p].id);E.x(other[0])
    else:E.cx(zc,index[p].id)
  for j,q in enumerate(index):
   if j!=p and code>>j&1:b.cx(index[p],q)
 ab.release_recycle('completion_index_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))
