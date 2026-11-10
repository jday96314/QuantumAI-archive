"""Low-workspace h=round(u^2/64), 0<=u<=256.

Split u=8v+r. Then h=v^2+floor((v*r+k(r))/4), where
k(r)=floor((r^2+32)/16) in {2,3,4,5}. Only two low carry
bits are needed; their measurement phase depends on five original bits.
The upper input bit is the isolated endpoint guard, not a signed bit.
"""
import json
from pathlib import Path
import numpy as np
import kickmix as km
from short_adder import add


def anf(values,n):
 c=list(map(int,values))
 for bit in range(n):
  for s in range(1<<n):
   if s>>bit&1:c[s]^=c[s^(1<<bit)]
 return c


def terms(b,src,target,coeffs,phase=False):
 # Share monomials separately in two small source halves. Cross-half
 # products are measured immediately; phase-only cross terms need CZ only.
 split=len(src)//2;needed=[set([0]),set([0])]
 for mon,value in enumerate(coeffs):
  if value:needed[0].add(mon&((1<<split)-1));needed[1].add(mon>>split)
 records=[];allocated=[]
 def bank(bits,masks):
  out={0:True}
  def get(mask):
   if mask in out:return out[mask]
   low=mask&-mask;j=low.bit_length()-1;rest=mask^low
   if not rest:out[mask]=bits[j]
   else:
    left=get(rest);q=b.alloc_qubits(1)[0];b.init_and(left,bits[j],q);records.append((left,bits[j],q));allocated.append(q);out[mask]=q
   return out[mask]
  for mask in sorted(masks):get(mask)
  return out
 lo=bank(src[:split],needed[0]);hi=bank(src[split:],needed[1])
 for mon,value in enumerate(coeffs):
  if not value:continue
  left=lo[mon&((1<<split)-1)];right=hi[mon>>split];tmp=None
  if not phase and left is not True and right is not True:
   tmp=b.alloc_qubits(1)[0];b.init_and(left,right,tmp)
  for j,q in enumerate(target):
   if not value>>j&1:continue
   if phase:b.push_condition(q);b.cz(left,right);b.pop_condition()
   elif tmp is not None:b.cx(tmp,q)
   elif left is True and right is True:b.x(q)
   else:b.cx(right if left is True else left,q)
  if tmp is not None:b.del_and(left,right,tmp);b.free(km.array([tmp]))
 for left,right,q in reversed(records):b.del_and(left,right,q);b.free(km.array([q]))


SQUARE=anf([v*v for v in range(32)],5)
LOW=anf([(((v>>3)*(v&7)+(((v&7)**2+32)//16))&3) for v in range(32)],5)


def apply(b,u,h,inverse=False):
 v=u[3:8];r=u[:3];crumb=b.alloc_qubits(2);acc=km.array(list(crumb)+list(h[:10]))
 def square():terms(b,v,h[:10],SQUARE);b.cx(u[8],h[10])
 def kinit(back=False):
  k=b.alloc_qubits(1)[0];b.init_and(r[1],r[2],k)
  def kbits():
   b.cx(r[2],crumb[0]);b.cx(k,crumb[0]);b.ccx(k,r[0],crumb[0]);b.x(crumb[1]);b.cx(k,crumb[1])
  if not back:kbits()
  add(b,km.array([k]),h[:10],subtract=back)
  if back:kbits()
  b.del_and(r[1],r[2],k);b.free(km.array([k]))
 def cross(back=False):
  for i in (reversed(range(3)) if back else range(3)):
   masked=b.alloc_qubits(5)
   for q,a in zip(masked,v):b.init_and(r[i],a,q)
   add(b,masked,acc[i:],subtract=back)
   for q,a in zip(masked,v):b.del_and(r[i],a,q)
   b.free(masked)
 if inverse:
  terms(b,u[:5],crumb,LOW);cross(True);kinit(True);square()
 else:
  square();kinit();cross();mask=b.alloc_bits(2);b.hmr(crumb,mask);terms(b,u[:5],mask,LOW,True);b.bit_store(mask,False);b.free(mask)
 b.free(crumb)


def build(inverse=False):
 b=km.CircuitBuilder();u=b.create_quantum_register(9,'u');h=b.create_quantum_register(11,'h');apply(b,u,h,inverse);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test():
 cf=build();ci=build(True);vals=np.resize(np.arange(257),1024);want=(vals*vals+32)//64;bits=lambda v,w:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool);rows=[]
 for seed in [142054,142055]:
  for rev,c in [(False,cf),(True,ci)]:
   sim=km.Simulator(batch_size=len(vals),seed=seed);sim.use_same_registers_as(c);regs=c.register_data;sim.write_across_shots(regs['u'],bits(vals,9))
   if rev:sim.write_across_shots(regs['h'],bits(want,11))
   sim.do(c);assert np.array_equal(sim.read_across_shots(regs['u']),bits(vals,9));bad=np.flatnonzero((sim.read_across_shots(regs['h'])!=bits(np.zeros_like(want) if rev else want,11)).any(0));phase=np.flatnonzero(sim.read_phase_flipped_across_shots());public=set(regs['u'])|set(regs['h']);sc=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(sc).any(0));r=dict(seed=seed,inverse=rev,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=bad[:8].tolist(),phase=phase[:8].tolist(),junk=junk[:8].tolist());rows.append(r);print(json.dumps(r),flush=True);assert not len(bad) and not len(phase) and not len(junk)
   sim.do(cf if rev else ci);assert np.array_equal(sim.read_across_shots(regs['u']),bits(vals,9));assert np.array_equal(sim.read_across_shots(regs['h']),bits(want if rev else np.zeros_like(want),11));assert not sim.read_phase_flipped_across_shots().any();assert not sim.read_across_shots(sc).any()
 for rev,c in [(False,cf),(True,ci)]:c.to_file(f'candidates/rounded_square_{"inverse" if rev else "forward"}.kmb',format='kmb')
 result=dict(status='PASS',shots=4*len(vals),rows=rows,scope='all 257 allowed magnitudes, endpoint 256, two seeds, fresh inverse and both compositions, clean scratch and phase');Path('results/rounded_square_native.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':test()
