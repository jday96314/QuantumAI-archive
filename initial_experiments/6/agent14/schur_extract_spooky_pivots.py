"""Quadratic Schur extraction with coherent adjacent-node pivot repair.

Consecutive swaps are suppressed, leaving a probability of order 1/q^2
for a non-completion scalar pivot. All node descriptors remain coherent.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import setup
from schur_extract import body,mapped,QRE,BRE,linear_map,imul_text,axpy_text
from schur_pivot_spooky import build as pivot
from spooky_extension_norm_inverse import build as inverse
from schur_control_templates import axpy,scale,flagger
from schur_shared_control import build as shared_control
from functools import lru_cache

@lru_cache(None)
def predicate(width,value):
 b=km.CircuitBuilder();index=b.create_quantum_register(width,'index');target=b.create_quantum_register(1,'target')[0];b.flip_if_equal(index,value,target=target);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))
from quadratic_schur_model import node_maps,transform,F,einv,bmul


def controls(N,i):
 n=N//2;out=[]
 if i:out.append((i-1,(2*i)^(2*(i-1))))
 if i<n-1:out.append((i,(2*i)^(2*(i+1))))
 return out


def build(N,reverse=False):
 start=time.time();n=N//2;m,k,fw,bk,tr=setup();nodes,kap,_,_=node_maps(n);cb=km.CircuitBuilder()
 for i in range(N):cb.create_quantum_register(12,f'S{i}')
 cb.create_quantum_register(12,'g');width=max(1,(n-1).bit_length());cb.create_quantum_register(width,'index');pairghostcount=24*n*(n-1)//2;ghostcount=pairghostcount+36*n;cb.create_classical_register(ghostcount,'ghost');headers=str(cb.finish_circuit())+'\n';data=[list(range(24*i,24*i+24)) for i in range(n)];g=list(range(12*N,12*N+12));index=list(range(12*N+12,12*N+12+width));work=list(range(index[-1]+1,index[-1]+4));shared=list(range(work[-1]+1,work[-1]+13));scratch=list(range(shared[-1]+1,shared[-1]+74));inv=inverse(reverse);ib=body(inv);ops=linear_map(N);pair_index={(i,j):sum(n-1-z for z in range(i))+j-i-1 for i in range(n) for j in range(i+1,n)};path=Path(f'candidates/schur_extract_spooky_pivots_n{N}_{"inverse" if reverse else "forward"}.kmx')
 with path.open('w') as f:
  f.write(headers)
  mask_live=False
  def emit(c,ports,pivot_index=None):
   size=len(ports);qm=dict(enumerate(ports));available=[q for q in (scratch if mask_live else shared+scratch) if q not in ports];qm.update({size+j:q for j,q in enumerate(available[:c.num_qubits-size])});f.write(BRE.sub(lambda x:f'b{pairghostcount+36*pivot_index+int(x[1]) if pivot_index is not None and int(x[1])<36 else ghostcount+int(x[1])-(36 if pivot_index is not None else 0)}',mapped(body(c),qm)))
  def swap(i):
   for a,b in zip(data[i],data[i+1]):f.write(f'CX q{a} q{b}\nCCX q{flags[i]} q{b} q{a}\nCX q{a} q{b}\n')
  def eq(value,target):emit(predicate(width,value),index+[target])
  def occupied():eq(0,work[2]);f.write(f'X q{work[2]}\n')
  def detect(i):
   prev=data[i-1][12:] if i else data[i][12:];emit(flagger(i==0),g+data[i][12:]+data[i+1][12:]+prev+[work[2],flags[i]])
  def store(i):
   for bit,q in enumerate(index):
    if (i+1)>>bit&1:f.write(f'CX q{flags[i]} q{q}\n')
  if not reverse:
   for op,a,b in ops:f.write(f'{op.upper()} q{a} q{b}\n')
  for i in (reversed(range(n)) if reverse else range(n)):
   cur=data[i];prev=data[i-1] if i else cur;cc=controls(N,i);flags={i-1:work[0],i:work[1]}
   if i:eq(i,flags[i-1])
   if reverse and i<n-1:eq(i+1,flags[i])
   def adjust_node():
    for flag,delta in cc:emit(axpy(tr(delta,fw)),[flags[flag]]+cur[12:]+cur[:12])
   def convert():emit(pivot(tr(nodes[i]&4095,fw),reverse,first=i==0),cur+g+prev,pivot_index=i)
   if reverse:convert();adjust_node()
   elif i<n-1:
    occupied();detect(i);store(i);occupied();f.write(f'CX q{flags[i]} q{work[2]}\n');swap(i)
   constants=tuple(tr(de,fw) for _,de in cc)
   shared_ports=[flags[fl] for fl,_ in cc]+cur[12:]+shared
   if cc and i<n-1:emit(shared_control(constants),shared_ports);mask_live=True
   for j in (reversed(range(i+1,n)) if reverse else range(i+1,n)):
    target=data[j];delta=tr((2*i)^(2*j),fw);constant=tr(kap[i]^kap[j],fw);active=[(fl,de) for fl,de in cc if not (fl==i and j==i+1)]
    def translate():
     f.write(''.join(f'CX q{a} q{b}\n' for a,b in zip(cur,target)));f.write(axpy_text(cur[12:],target[:12],delta,m))
     if cc:
      f.write(''.join(f'CX q{a} q{b}\n' for a,b in zip(shared,target[:12])))
      if j==i+1:emit(axpy(tr(cc[-1][1],fw)),[flags[i]]+cur[12:]+target[:12])
    def multiply(back=False):
     if not back:
      for word in [target[:12],target[12:]]:f.write(imul_text(word,constant,m))
     for flag,de in (reversed(active) if back else active):
      dd=tr(F.mul(de,de)^de,fw);rho=m.gmul(constant^dd,m.ginv(constant))
      for word in [target[:12],target[12:]]:emit(scale(rho,back),[flags[flag]]+word)
     if back:
      for word in [target[:12],target[12:]]:f.write(imul_text(word,constant,m,True))
    if reverse:multiply(True)
    else:translate()
    qm=dict(enumerate(target));qm.update({24+x:scratch[x] for x in range(inv.num_qubits-24)});txt=mapped(ib,qm);base=24*pair_index[i,j];txt=BRE.sub(lambda x:f'b{base+int(x[1]) if int(x[1])<24 else ghostcount+int(x[1])-24}',txt);f.write(txt)
    if reverse:translate()
    else:multiply()
   if cc and i<n-1:emit(shared_control(constants,True),shared_ports);mask_live=False
   if reverse:
    if i<n-1:swap(i);store(i);occupied();detect(i);occupied()
   else:
    adjust_node();convert()
    if i<n-1:eq(i+1,flags[i])
   if i:eq(i,flags[i-1])
  if reverse:
   for op,a,b in reversed(ops):f.write(f'{op.upper()} q{a} q{b}\n')
 c=km.Circuit.from_file(str(path),format='kmx');c.to_file(str(path.with_suffix('.kmb')),format='kmb');r=dict(N=N,reverse=reverse,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,seconds=time.time()-start,scope='one adjacent swap recorded by a coherent index; later scalar pivots unsupported; retained X outcomes until reverse');Path(f'results/schur_extract_spooky_pivots_n{N}_{"inverse" if reverse else "forward"}_build.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True);return c


def model(s):
 n=len(s)//2;nodes,kap,P,maps=node_maps(n);v=transform(s,maps);nodes=list(nodes);kap=list(kap);g=1;A=B=0;out=[];flags=[];ghost=[];pivotghost=[]
 for i in range(n):
  curbf=F.mul(g,v[i]>>12)^B
  if i<n-1:
   nxtbf=F.mul(g,v[i+1]>>12)^B;ff=not any(flags) and not curbf and bool(nxtbf);flags.append(int(ff))
   if ff:v[i],v[i+1]=v[i+1],v[i];nodes[i],nodes[i+1]=nodes[i+1],nodes[i];kap[i],kap[i+1]=kap[i+1],kap[i]
  b=v[i]>>12;a=(v[i]&4095)^F.mul(nodes[i]&4095,b);pivotghost.append(a|(b<<12)|(g<<24))
  for j in range(i+1,n):
   den=v[j]^a^bmul(nodes[j],b);ghost.append(den);v[j]=bmul(einv(den),kap[i]^kap[j]) if den else 0
  aa=F.mul(g,a)^A;bb=F.mul(g,b)^B;out.extend([aa,bb]);g=F.mul(F.mul(bb or 1,bb or 1),F.inv(g));A=aa^bb;B=bb
 return out,g,flags,ghost,pivotghost


def test(N=6,shots=512):
 cf=build(N);ci=build(N,True);m,k,fw,bk,tr=setup();rng=np.random.default_rng(142052+N);ss=rng.integers(0,4096,(shots,N));n=N//2
 for j in range(0,shots,2):
  d=(j//2)%(n+1);a=rng.choice(np.arange(1,4096),d,replace=False).tolist();e=rng.integers(1,4096,d).tolist();ss[j]=F.syndrome(a,e,N)
 # Force swaps in many arbitrary-data cases by making the first raw high
 # sample zero and the next nonzero, using the inverse binary transform.
 from tower_inverse import map_ops
 # Random arbitrary syndromes already exercise full basis permutations;
 # an additional forced distribution is prepared through the host map below.
 _,_,_,maps=node_maps(n);cols=[]
 for j in range(N):
  for bit in range(12):
   x=1<<bit;col=0
   for ii,qq in enumerate(maps):col^=bmul(qq[j],x)<<(24*ii)
   cols.append(col)
 ops=map_ops(cols)
 for case in range(1,shots,4):
  vals=[int(x) for x in rng.integers(0,1<<24,n)];vals[0]&=4095;vals[1]|=1<<12;packed=sum(x<<(24*i) for i,x in enumerate(vals))
  for op,a,b in reversed(ops):
   if op=='swap':
    if ((packed>>a)^(packed>>b))&1:packed^=(1<<a)|(1<<b)
   elif packed>>a&1:packed^=1<<b
  ss[case]=[(packed>>(12*i))&4095 for i in range(N)]
 modelrows=[model(s.tolist()) for s in ss];want=np.array([[tr(x,fw) for x in v[0]] for v in modelrows]);gv=np.array([tr(v[1],fw) for v in modelrows]);flagrows=np.array([v[2] for v in modelrows],dtype=bool).T;indices=np.array([next((i+1 for i,f in enumerate(v[2]) if f),0) for v in modelrows]);width=max(1,(n-1).bit_length());ff=((indices[None,:]>>np.arange(width)[:,None])&1).astype(bool);gg=np.array([[tr(x&4095,fw)|(tr(x>>12,fw)<<12) for x in v[3]] for v in modelrows]).T;bits=lambda v:((v[None,:]>>np.arange(12)[:,None])&1).astype(bool);reports=[]
 for reverse,c,src,dst,gin,gout in [(False,cf,ss,want,np.ones(shots,dtype=int),gv),(True,ci,want,ss,gv,np.ones(shots,dtype=int))]:
  sim=km.Simulator(batch_size=shots,seed=142052+reverse);sim.use_same_registers_as(c)
  for i in range(N):sim.write_across_shots(c.register_data[f'S{i}'],bits(src[:,i]))
  sim.write_across_shots(c.register_data['g'],bits(gin))
  if reverse:sim.write_across_shots(c.register_data['index'],ff)
  sim.do(c);errors={}
  for i in range(N):
   bad=np.flatnonzero((sim.read_across_shots(c.register_data[f'S{i}'])!=bits(dst[:,i])).any(0))
   if len(bad):errors[f'S{i}']=bad[:8].tolist()
  bad=np.flatnonzero((sim.read_across_shots(c.register_data['g'])!=bits(gout)).any(0))
  if len(bad):errors['g']=bad[:8].tolist()
  assert np.array_equal(sim.read_across_shots(c.register_data['index']),np.zeros_like(ff) if reverse else ff)
  public={km.q(i) for i in range(12*N+12+width)};sc=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);junk=np.flatnonzero(sim.read_across_shots(sc).any(0))
  all_masks=sim.read_across_shots(c.register_data['ghost']);pairghostcount=24*n*(n-1)//2;masks=all_masks[:pairghostcount].reshape(-1,24,shots);gg=np.array([[m.gmul(int(x)&4095,int(x)&4095)^m.gmul(int(x)&4095,int(x)>>12)^m.gmul(k,m.gmul(int(x)>>12,int(x)>>12)) for x in row] for row in gg]);gbits=((gg[:,None,:]>>np.arange(24)[None,:,None])&1).astype(bool);expected=np.zeros(shots,dtype=bool) if reverse else ((masks&gbits).sum((0,1))%2).astype(bool);pv=np.array([[tr((x>>(12*j))&4095,fw) for x in v[4] for j in range(3)] for v in modelrows]).reshape(shots,n,3);pbits=((pv.transpose(1,2,0)[:,:,None,:]>>np.arange(12)[None,None,:,None])&1).astype(bool).reshape(n,36,shots);expected^=np.zeros(shots,dtype=bool) if reverse else ((all_masks[pairghostcount:].reshape(n,36,shots)&pbits).sum((0,1))%2).astype(bool);phase=np.flatnonzero(sim.read_phase_flipped_across_shots()!=expected);r=dict(reverse=reverse,errors=errors,junk=junk[:8].tolist(),phase=phase[:8].tolist());reports.append(r);print(json.dumps(r),flush=True);assert not errors and not len(junk) and not len(phase)
  if not reverse:
   sim.do(ci)
   for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(src[:,i]))
   assert np.array_equal(sim.read_across_shots(c.register_data['g']),bits(gin));assert not sim.read_across_shots(c.register_data['index']).any();assert not sim.read_across_shots(c.register_data['ghost']).any();assert not sim.read_across_shots(sc).any();assert not sim.read_phase_flipped_across_shots().any()
 result=dict(status='PASS',N=N,shots=shots,swaps=int(flagrows.sum()),rows=reports,scope='forced first swap, actual syndromes and arbitrary inputs; exact pivot/frame/index, forward ghost phase, inverse and pair cleanup; only one repair allowed');Path(f'results/schur_extract_spooky_pivots_n{N}_native.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);p.add_argument('--shots',type=int,default=512);a=p.parse_args();test(a.n,a.shots)
