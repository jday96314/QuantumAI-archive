"""Two-payload inverse Bezout transform through a generic Euclidean chart.

Sources are restored. Zero pivots are reversibly armed, defining a supported
permutation but not the desired Bezout map on degree-gap charts.
"""
from contextlib import contextmanager
import argparse,json,time
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from low_workspace_inverse import bind
from fast_dqi import FastBuilder


def divide(m,E,f,g,inverse=False,one=False,provided_inverse=None):
 n=len(g)-1;mm=len(f)-1
 @contextmanager
 def inv():
  if one:yield None
  elif provided_inverse is not None:yield provided_inverse
  else:
   with m.inv_word(E,g[-1]) as gi:yield gi
 with inv() as gi,m.opened(E,'bez_source_divide',12) as tmp:
  for i in (range(n,mm+1) if inverse else reversed(range(n,mm+1))):
   if not inverse and not one:m.mon_slot(E,f[i],gi,g[-1],tmp)
   for j in range(n):E.mac(m.prod_rt(f[i-n+j],f[i],g[j]))
   if inverse and not one:m.mon_slot(E,f[i],g[-1],gi,tmp)


def apply(m,E,f,g,payloads,inverse=False,g_one=False):
 mm,n=len(f)-1,len(g)-1;assert all(len(v)==mm for v in payloads) and mm>=n>=0
 if mm==n:
  divide(m,E,f,g,one=g_one);apply(m,E,g,f[:n],payloads,inverse);divide(m,E,f,g,True,g_one);return
 assert mm==n+1
 with m.opened(E,'bez_zero',1) as zz:
  zero=zz[0]
  with m.eq_test(E,g[-1],0) as z:E.cx(z,zero)
  E.cx(zero,g[-1][0])
  def scale_last(scale,eraser):
   with m.opened(E,'bez_fast_scale',12) as tmp:
    for y in payloads:m.mon_slot(E,y[-1],scale,eraser,tmp)
  def cross(post=False):
   for i in range(n):
    for y in payloads:E.mac(m.prod_rt(y[-1] if post else y[i],g[i],y[i] if post else y[-1]))
  if not n:
   with m.opened(E,'bez_fast_pivot',12) as d:
    rt=m.prod_rt(d,f[-1],g[-1]);E.mac(rt)
    with m.inv_word(E,d) as di:scale_last(d if inverse else di,di if inverse else d)
    E.free(rt)
  else:
   with m.inv_word(E,g[-1]) as gi:
    scale_last(g[-1] if inverse else gi,gi if inverse else g[-1]);cross(inverse);divide(m,E,f,g,provided_inverse=gi)
   if inverse:apply(m,E,g,f[:n],[y[:n] for y in payloads],True)
   with m.inv_word(E,f[-1]) as qi:scale_last(f[-1] if inverse else qi,qi if inverse else f[-1])
   if not inverse:apply(m,E,g,f[:n],[y[:n] for y in payloads])
   with m.inv_word(E,g[-1]) as gi:
    divide(m,E,f,g,True,provided_inverse=gi);cross(not inverse);scale_last(g[-1] if inverse else gi,gi if inverse else g[-1])
  E.cx(zero,g[-1][0])
  with m.eq_test(E,g[-1],0) as z:E.cx(z,zero)


def build(t,inverse=False):
 m=setup()[0];m.inv_word=bind(m);b=FastBuilder(km.CircuitBuilder());ff=[b.create_quantum_register(12,f'f{i}') for i in range(t+1)];gg=[b.create_quantum_register(12,f'g{i}') for i in range(t)];ys=[[b.create_quantum_register(12,f'y{p}_{i}') for i in range(t)] for p in range(2)];ab=Adapter(b);sc=ab.qubits(3);E=m.Em(ab,*sc,m.Alloc(ab));lead=E.A.open('bez_monic_lead',12);E.x(lead[0]);apply(m,E,[[q.id for q in w] for w in ff],[[q.id for q in w] for w in gg]+[lead],[[[q.id for q in w] for w in bank] for bank in ys],inverse,True);E.x(lead[0]);E.A.close('bez_monic_lead',lead);ab.release_recycle('bez_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def test(t=3,shots=512):
 from derivative_inverse_model import bezout,solve,matvec,divide as hostdiv
 m,_,fw,bk,tr=setup();from quadratic_schur_model import F
 rng=np.random.default_rng(142085+t);samples=[];rejected=0
 while len(samples)<shots:
  f=rng.integers(0,4096,t+1).tolist();g=rng.integers(0,4096,t).tolist()+[1];a,b=f,g;good=True
  while len(b)>1:
   _,r=hostdiv(a,b)
   if len(r)!=len(b)-1:good=False;break
   a,b=b,r
  if not good or not b or not b[0]:rejected+=1;continue
  y=rng.integers(0,4096,(2,t)).tolist();Z=bezout(f,g);want=[solve(Z,v) for v in y];samples.append((f,g[:-1],y,want))
 bits=lambda v:((np.asarray(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool);cf=build(t);ci=build(t,True);rows=[]
 for inverse,c in [(False,cf),(True,ci)]:
  sim=km.Simulator(batch_size=shots,seed=142085+inverse);sim.use_same_registers_as(c)
  for key,field in [('f',0),('g',1)]:
   for i in range(t+int(key=='f')):sim.write_across_shots(c.register_data[f'{key}{i}'],bits([tr(v[field][i],fw) for v in samples]))
  for p in range(2):
   for i in range(t):sim.write_across_shots(c.register_data[f'y{p}_{i}'],bits([tr(v[3 if inverse else 2][p][i],fw) for v in samples]))
  sim.do(c);errors={}
  for key,field in [('f',0),('g',1)]:
   for i in range(t+int(key=='f')):
    if not np.array_equal(sim.read_across_shots(c.register_data[f'{key}{i}']),bits([tr(v[field][i],fw) for v in samples])):errors[f'{key}{i}']=True
  for p in range(2):
   for i in range(t):
    if not np.array_equal(sim.read_across_shots(c.register_data[f'y{p}_{i}']),bits([tr(v[2 if inverse else 3][p][i],fw) for v in samples])):errors[f'y{p}_{i}']=True
  public={q for reg in c.register_data.values() for q in reg};junk=km.array([km.q(i) for i in range(c.num_qubits) if km.q(i) not in public]);r=dict(inverse=inverse,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,scratch=bool(sim.read_across_shots(junk).any()),phase=bool(sim.read_phase_flipped_across_shots().any()));rows.append(r);print(json.dumps(r),flush=True);assert not errors and not r['scratch'] and not r['phase'];sim.do(cf if inverse else ci);assert not sim.read_across_shots(junk).any();assert not sim.read_phase_flipped_across_shots().any()
 result=dict(status='PASS',t=t,shots=shots*2,rows=rows,rejected_generic_charts=rejected,scope='two arbitrary payloads, generic Euclidean degree chart, inverse Bezout matrix action and literal inverse, restored sources/scratch/phase');Path(f'results/bezout_payload_fast_t{t}_native.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--t',type=int,default=3);p.add_argument('--shots',type=int,default=256);a=p.parse_args();test(a.t,a.shots)
