"""Root-gated native scan for a monic proper quadratic-Schur rational pair."""
from functools import lru_cache
from types import SimpleNamespace
from pathlib import Path
import json
import numpy as np
import kickmix as km
from fast_dqi import FastBuilder
from recursive_forney import RecursiveForney
from extension_inverse import setup
from schur_synthesize import plan
from schur_extract import body,mapped,BRE,imul_text,axpy_text
from tower_inverse import map_ops
from quadratic_schur_model import node_maps,F

@lru_cache(None)
def phase_template():
 b=FastBuilder(km.CircuitBuilder());s=b.create_quantum_register(12,'s');d=b.create_quantum_register(12,'d');o=b.create_quantum_register(12,'o');mask=b.create_classical_register(12,'mask')
 class Kernel(RecursiveForney):gf=km.GF2Field(12)
 Kernel()._forney_kernel(b,s,d,o,mask);c=b.finish_circuit()
 return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(12,c.num_bits)))


def write_scan(f,N,M,private_base,repair=False,banks=None,scratch_base=None,private_qubits=None):
 m,k,fw,bk,tr=setup();steps,nn,bb,free=plan(N) if not repair else __import__('schur_synthesize_repair').plan(N)
 if banks is not None:nn,bb=banks
 n=N//2;extra_flags=n-1 if repair else 0;gf=km.GF2Field(12);_,_,P,_=node_maps(n);host=gf;linearfield=SimpleNamespace(gmul=gf.mul);templ=phase_template();lines=body(templ)
 base=12*N+36+extra_flags if scratch_base is None else scratch_base
 scratch=list(range(base,base+templ.num_qubits-36)) if private_qubits is None else list(private_qubits[:templ.num_qubits-36])
 assert len(scratch)==templ.num_qubits-36 and not (set(scratch)&{q for w in nn+bb for q in w})
 qmap={**dict(enumerate(bb[0]+bb[1]+nn[0])),**{36+j:q for j,q in enumerate(scratch)}};kernel=mapped(lines,qmap)
 def cmap(word,ops):
  for op,a,b in ops:f.write(f'{op.upper()} q{word[a]} q{word[b]}\n')
 for word in nn+bb:cmap(word,map_ops(bk))
 odd=bb[3::2];even=bb[2::2];others=nn[1:];a=1;alpha=int(gf.primitive_element)
 for j in range(M):
  top=int(gf.pow(a,n));p=F.evaluate(P,a);assert p
  # Derivative times a includes the implicit leading term, since n is odd.
  # For even n it contributes only to the root accumulator.
  def evals(back=False):
   ops=[]
   for word in odd:ops.extend(f'CX q{x} q{y}\n' for x,y in zip(word,bb[1]))
   if n%2:
    ops.extend(f'X q{bb[1][bit]}\n' for bit in range(12) if top>>bit&1)
   for word in [bb[1]]+even:ops.extend(f'CX q{x} q{y}\n' for x,y in zip(word,bb[0]))
   if not n%2:ops.extend(f'X q{bb[0][bit]}\n' for bit in range(12) if top>>bit&1)
   for word in others:ops.extend(f'CX q{x} q{y}\n' for x,y in zip(word,nn[0]))
   f.write(''.join(reversed(ops) if back else ops))
  evals();f.write(f'INFO_TAG "EVAL_{j}"\n')
  invp=gf.pow(p,4094);cols=[gf.mul(invp,1<<i) for i in range(12)];trans=tuple(sum(((cols[i]>>b)&1)<<i for i in range(12)) for b in range(12));maskops=map_ops(trans)
  for op,x,y in maskops:
   for x,y in ([(x,y)] if op=='cx' else [(x,y),(y,x),(x,y)]):f.write(f'BIT_INVERT b{12*j+y} if b{12*j+x}\n')
  f.write(BRE.sub(lambda x:f'b{12*j+int(x[1]) if int(x[1])<12 else private_base+int(x[1])-12}',kernel))
  for op,x,y in reversed(maskops):
   for x,y in ([(x,y)] if op=='cx' else [(x,y),(y,x),(x,y)]):f.write(f'BIT_INVERT b{12*j+y} if b{12*j+x}\n')
  evals(True)
  if j+1<M:
   for bank in [nn,bb]:
    for d,word in enumerate(bank[1:],1):f.write(imul_text(word,gf.pow(alpha,d),linearfield))
   a=int(gf.mul(a,alpha))
 if M>1:
  for bank in [nn,bb]:
   for d,word in enumerate(bank[1:],1):f.write(imul_text(word,gf.pow(alpha,d*(M-1)%4095),linearfield,True))
 for word in nn+bb:cmap(word,reversed(map_ops(bk)))
 return dict(kernel_Q=templ.num_qubits,kernel_T=templ.operation_statistics().num_toffoli_like)


def test():
 from io import StringIO
 N=6;M=63;cb=km.CircuitBuilder()
 for i in range(N):cb.create_quantum_register(12,f'S{i}')
 cb.create_quantum_register(12,'g')
 for j in range(M):cb.create_classical_register(12,f'E{j}')
 headers=str(cb.finish_circuit());f=StringIO();f.write(headers+'\n');r=write_scan(f,N,M,12*M);c=km.Circuit(f.getvalue());m,k,fw,bk,tr=setup();steps,nn,bb,free=plan(N);rng=np.random.default_rng(142048);shots=512;_,_,P,_=node_maps(N//2);mask=rng.integers(0,4096,(shots,M));pos=[];val=[];nums=[];dens=[]
 for s in range(shots):
  d=s%(N//2+1);aa=rng.choice(np.arange(M),d,replace=False).tolist();xs=[F.ex[j] for j in aa];es=rng.integers(1,4096,d).tolist();B=[1]
  for a in xs:B=F.product(B,[a,1])
  Np=[0]
  for a,e in zip(xs,es):
   q=[0]*(len(B)-1);q[-1]=B[-1]
   for i in reversed(range(len(q)-1)):q[i]=B[i+1]^F.mul(a,q[i+1])
   Np=F.add(Np,F.scale(q,F.mul(e,F.mul(a,F.evaluate(P,a)))))
  for a,e in zip(xs,es):
   got=F.mul(F.evaluate(Np,a),F.inv(F.mul(F.evaluate(F.derivative(B),a),F.mul(a,F.evaluate(P,a)))))
   assert got==e,(a,e,got,Np,B)
  pad=N//2-d;Bp=[0]*pad+B;Np=[0]*pad+Np;Bp=(Bp+[0]*(N//2+1))[:N//2];Np=(Np+[0]*(N//2))[:N//2];nums.append([tr(x,fw) for x in Np]);dens.append([tr(x,fw) for x in Bp]);pos.append(aa);val.append(es)
 nums=np.array(nums);dens=np.array(dens);bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);sim=km.Simulator(batch_size=shots,seed=142048);sim.use_same_registers_as(c)
 for ws,vs in [(nn,nums),(bb,dens)]:
  for j,w in enumerate(ws):sim.write_across_shots(km.array([km.q(q) for q in w]),bits(vs[:,j]))
 for j in range(M):sim.write_across_shots(c.register_data[f'E{j}'],bits(mask[:,j]))
 sim.do(c);expected=np.array([sum((int(mask[i,j])&e).bit_count() for j,e in zip(pos[i],val[i]))%2 for i in range(shots)],dtype=bool);bad=np.flatnonzero(sim.read_phase_flipped_across_shots()!=expected);print('phase',bad[:10].tolist(),flush=True)
 assert not len(bad)
 for ws,vs in [(nn,nums),(bb,dens)]:
  for j,w in enumerate(ws):assert np.array_equal(sim.read_across_shots(km.array([km.q(q) for q in w])),bits(vs[:,j]))
 for j in range(M):assert np.array_equal(sim.read_across_shots(c.register_data[f'E{j}']),bits(mask[:,j]))
 used={q for w in nn+bb for q in w};sc=km.array([km.q(q) for q in range(c.num_qubits) if q not in used]);assert not sim.read_across_shots(sc).any();r.update(status='PASS',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=shots,scope='constructed error patterns, all supported weights, implicit monic polynomial padding, restored coefficients/masks/scratch, exact error-magnitude phase');Path('results/schur_scan_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':test()
