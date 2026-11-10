"""Complete two-stage fused inverse-square-root prototype.

Uses the exhaustively screened extra=1 independently rounded Booth rows.
This is a resource/correctness experiment until a complete verifier accepts it.
"""
import argparse,json,math
from fractions import Fraction
from pathlib import Path
import numpy as np
import kickmix as km
from binary_centered_invsqrt import cshift as input_shift,del_flag_if_zero
from short_adder import add
from booth_adder import accumulate
from rounded_coordinate_inplace import apply as square_coordinate
from streamed_unlookup import erase
from invsqrt.nr_reference_circuit import _phase_by_bits

TABLE=np.load('parameters/fused_q5_e1_booth.npy')
A=(961<<17,1217<<17)


def load(b,address,regs,values):
 def leaf(row,control):
  for r,v in zip(regs,row):b.gf2_iadd(int(v),target=r,control=control)
 def visit(control,level,lo,hi):
  if np.all(values[lo:hi]==values[lo]):leaf(values[lo],control);return
  bit=address[level];mid=(lo+hi)//2;q=b.alloc_qubits(1)[0];b.x(bit);b.init_and(control,bit,q);b.x(bit);visit(q,level-1,lo,mid);b.cx(control,q);visit(q,level-1,mid,hi);b.del_and(control,bit,q);b.free(km.array([q]))
 mid=len(values)//2;b.x(address[-1]);visit(address[-1],len(address)-2,0,mid);b.x(address[-1]);visit(address[-1],len(address)-2,mid,len(values))


def shift(b,y,flags,inverse=False):
 occ=[28,29,31,35]
 for k in (reversed(range(4)) if inverse else range(4)):
  step=1<<k;indices=[i for i in range(min(occ[k],len(y))) if i+step<len(y)]
  for i in (indices if inverse else reversed(indices)):
   a,c=y[i],y[i+step]
   if i==27 and k==0:b.cx(flags[k+1],a);b.cx(flags[k+1],c)
   elif i+step>=occ[k]:
    if inverse:b.cx(c,a);b.del_and(flags[k+1],a,c)
    else:b.init_and(flags[k+1],a,c);b.cx(c,a)
   else:b.cswap(flags[k+1],a,c)


def build(limit=90,product=14):
 b=km.CircuitBuilder();x=b.create_quantum_register(22,'x');flags=b.alloc_qubits(5)
 def normalize(inverse=False):
  for k in (range(5) if inverse else reversed(range(5))):
   s=1<<k
   if inverse:input_shift(b,x,s,flags[k],22-s,True);del_flag_if_zero(b,x[-s:],flags[k])
   else:b.flip_if_equal(x[-s:],0,target=flags[k]);input_shift(b,x,s,flags[k],22-s)
 b.info_tag('normalize_forward');normalize();b.info_tag('square_forward');h=b.alloc_qubits(11);square_coordinate(b,x,h);y0=b.alloc_qubits(24);curve=b.alloc_qubits(7);address=km.array(list(x[13:21])+[flags[0]]);idx=np.arange(512);coarse=(idx&255)>>1;parity=idx>>8;rem=(TABLE[:,0]+(1<<25)+(coarse<<18)+(1<<17)-np.array(A)[parity])<<1;assert rem.min()>=0 and rem.max()<(1<<24);b.info_tag('coefficient_lookup');load(b,address,[y0,curve],np.stack([rem,TABLE[:,2]],axis=1));extra=b.alloc_qubits(3);y=km.array(list(y0)+list(extra)+[x[21]]);cc=b.alloc_qubits(1);lc=None
 def high_integrate(c,inverse=False,signed=False):
  high=y[19:27]
  if signed:b.cx(x[13],high)
  add(b,c,high,subtract=inverse,qubit_limit=limit)
  if signed:b.cx(x[13],high)
 def curvature(inverse=False):
  if inverse:high_integrate(cc,True)
  b.x(x[21]);accumulate(b,curve,km.array(list(h)+[x[21]]),km.array(list(y[:19])+list(cc)),-6,subtract=inverse,rounded=True,borrow=y[26:27],auxiliary=y[24:26],qubit_limit=limit,product_limit=product);b.x(x[21])
  if not inverse:high_integrate(cc)
 b.info_tag('curve_forward');curvature();b.info_tag('release_curve');hghost=b.alloc_bits(11);cghost=b.alloc_bits(7);b.hmr(h,hghost);b.hmr(curve,cghost);b.free(h);b.free(curve);slope=b.alloc_qubits(17);b.info_tag('slope_lookup');load(b,address,[slope],TABLE[:,1:2]);lc=b.alloc_qubits(1)
 def linear(inverse=False):
  if inverse:high_integrate(lc,True,True)
  b.x(x[13]);b.x(x[21]);accumulate(b,slope,x[:14],km.array(list(y[:19])+list(lc)),-11,subtract=not inverse,rounded=True,borrow=y[24:27],auxiliary=x[21:22],qubit_limit=limit,product_limit=product);b.x(x[21]);b.x(x[13])
  if not inverse:high_integrate(lc,False,True)
 def affine(inverse=False):
  def constant(back):
   a,c,d=y[24:27];q=b.alloc_qubits(1)[0]
   if not back:b.x(a);b.cx(a,c)
   b.init_and(a,c,q);b.x(q);b.cx(q,d);b.x(q);b.del_and(a,c,q)
   if back:b.cx(a,c);b.x(a)
   b.free(km.array([q]))
  def par(back):b.cx(flags[0],y[26])
  def coarse_shift(back):add(b,x[14:21],y[19:27],subtract=not back,qubit_limit=limit)
  for f in ([coarse_shift,par,constant] if inverse else [constant,par,coarse_shift]):f(inverse)
 b.info_tag('linear_forward');linear();b.info_tag('affine_forward');affine();b.info_tag('phase');head=b.alloc_qubits(10);phase_y=km.array(list(y)+list(head));shift(b,phase_y,flags);_phase_by_bits(b,phase_y);shift(b,phase_y,flags,True);b.free(head);b.info_tag('affine_reverse');affine(True);b.info_tag('linear_reverse');linear(True);b.free(lc);b.info_tag('curve_replay');sghost=b.alloc_bits(17);b.hmr(slope,sghost);b.free(slope);curve=b.alloc_qubits(7);load(b,address,[curve],TABLE[:,2:3]);h=b.alloc_qubits(11);square_coordinate(b,x,h)
 for q,c in zip(curve,cghost):b.push_condition(c);b.z(q);b.pop_condition()
 for q,c in zip(h,hghost):b.push_condition(c);b.z(q);b.pop_condition()
 b.bit_store(cghost,False);b.bit_store(hghost,False);b.free(cghost);b.free(hghost);b.info_tag('curve_reverse');curvature(True);b.free(cc);b.info_tag('square_erase');square_coordinate(b,x,h,True);b.free(h);packed=[int(r)|(int(c)<<24)|(int(s)<<31) for r,s,c in zip(rem,TABLE[:,1],TABLE[:,2])];b.info_tag('coefficient_erase');erase(b,address,km.array(list(y0)+list(curve)),packed,extra_outcomes=sghost,split=4);b.free(y0);b.free(curve);b.free(extra);b.info_tag('normalize_reverse');normalize(True);b.free(flags);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def model(xs):
 xs=np.asarray(xs,dtype=np.int64);k=21-np.floor(np.log2(xs)).astype(np.int64);mant=xs<<k;s=k//2;address=((mant>>13)&255)+((k&1)<<8);c0,c1,c2=TABLE[address].T;z=(mant&16383)^8192;signed=(mant&16383)-8192;u=(signed+16)>>5;h=(u*u+32)>>6
 def booth(a,z,L,e):
  out=np.zeros_like(z)
  for g in range((L+1)//2):
   bit=lambda i:0 if i<0 else (z>>min(i,L-1))&1
   digit=bit(2*g-1)+bit(2*g)-2*bit(2*g+1)
   for mag in (1,2):
    ee=e+2*g+int(mag==2);v=a<<ee if ee>=0 else (a>>-ee)+((a>>(-ee-1))&1);out+=np.where(abs(digit)==mag,np.sign(digit)*v,0)
  return out
 p=(c0<<1)-booth(c1,z,14,-11)+booth(c2,h,12,-6);return p.astype(np.uint64)<<(s+26).astype(np.uint64)


def test(c,shots=8192):
 rng=np.random.default_rng(142082);xs=rng.integers(1,3133453,shots);xs[:100]=np.r_[np.arange(1,65),np.array([max(1,min(3133452,(1<<i)+j)) for i in range(10,22) for j in (-1,0,1)])];sim=km.Simulator(batch_size=shots,seed=142082);sim.use_same_registers_as(c);reg=c.register_data['x'];bits=lambda v,w:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool);sim.write_across_shots(reg,bits(xs,22));sim.do(c);got=np.array([int(sim.read_shot_phase(i)*(1<<63))&((1<<64)-1) for i in range(shots)],dtype=np.uint64);want=model(xs);bad=np.flatnonzero(got!=want);junk=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in set(reg)]);r=dict(status='PASS' if not len(bad) and np.array_equal(sim.read_across_shots(reg),bits(xs,22)) and not sim.read_across_shots(junk).any() else 'FAIL',shots=shots,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,phase_errors=bad[:8].tolist(),scratch_qubits=np.flatnonzero(sim.read_across_shots(junk).any(1))[:8].tolist());print(json.dumps(r),flush=True);Path('results/fused_invsqrt_native.json').write_text(json.dumps(r,indent=2)+'\n');assert r['status']=='PASS'

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=90);p.add_argument('--product',type=int,default=14);a=p.parse_args();c=build(a.limit,a.product);c.to_file('candidates/fused_invsqrt_prototype.kmb',format='kmb');print(json.dumps(dict(Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like)),flush=True);test(c)
