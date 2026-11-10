"""A normalized projective jet endpoint with no pending measurement phases."""
import argparse,json,time
from pathlib import Path
import numpy as np
import kickmix as km
from projective_jet_step import build as jet_step
from schur_extract import body,mapped,BRE
from raw_cf_endpoint import test as old_test
from raw_cf_model import extract,endpoint
from quadratic_schur_model import F
from extension_inverse import setup


def build(N,inverse=False):
 n=N//2;b=km.CircuitBuilder()
 for i in range(N):b.create_quantum_register(12,f'S{i}')
 b.create_quantum_register(49,'endpoint_ghost');tail=[list(range(12*N+12*i,12*N+12*i+12)) for i in range(3)];tag=12*N+36
 # Text register declarations do not reserve the zero tail workspace.
 headers=str(b.finish_circuit())+'\n'+''.join(''.join(f'APPEND_TO_REGISTER q{q} r{N+1+i}\n' for q in w)+f'REGISTER r{N+1+i} "{name}"\n' for i,(name,w) in enumerate(zip(['tailN0','tailN1','tailB0'],tail)))
 c=jet_step(inverse);f=[headers];pool=tail[2]+list(range(tag+1,tag+80))
 def add_gate(op,*qs):f.append(op+' '+ ' '.join(f'q{q}' for q in qs)+'\n')
 def boundary():
  if inverse:add_gate('X',tag);add_gate('CX',tail[2][0],tag);add_gate('X',tail[2][0]);add_gate('CX',tag,tail[2][0]);add_gate('CX',tag,tail[0][0])
  else:add_gate('CX',tag,tail[0][0]);add_gate('X',tail[2][0]);add_gate('CX',tag,tail[2][0]);add_gate('CX',tail[2][0],tag);add_gate('X',tag)
 if inverse:boundary()
 for i in (reversed(range(n)) if inverse else range(n)):
  ports={'a':list(range(24*i,24*i+12)),'b':list(range(24*i+12,24*i+24)),'h':tail[1],'k':tail[0],'infinity':[tag]};qm={q.id:v for name,w in ports.items() for q,v in zip(c.register_data[name],w)};missing=[q for q in range(c.num_qubits) if q not in qm];qm.update(dict(zip(missing,pool)));f.append(BRE.sub(lambda x:f'b{49+int(x[1])}',mapped(body(c),qm)))
 if not inverse:boundary()
 return km.Circuit(''.join(f))


def test(N=6,shots=512):
 m,_,fw,bk,tr=setup();n=N//2;rng=np.random.default_rng(142119+N);src=[];want=[]
 for case in range(shots):
  d=case%(n+1);pos=rng.choice(np.arange(1,4096),d,replace=False).tolist();val=rng.integers(1,4096,d).tolist();hh,_=extract(F.syndrome(pos,val,N));a,b,e,j=endpoint(hh)
  if e:continue
  a=a+[0]*(2-len(a));src.append([tr(x,fw) for h in hh for x in h]);want.append([tr(a[0],fw),tr(a[1],fw),int(bool(b))])
 src=np.array(src);want=np.array(want);shots=len(src);bits=lambda v:((np.asarray(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool);cf=build(N);ci=build(N,True);rows=[]
 for inv,c in [(False,cf),(True,ci)]:
  sim=km.Simulator(batch_size=shots,seed=142119+inv);sim.use_same_registers_as(c)
  for i in range(N):sim.write_across_shots(c.register_data[f'S{i}'],bits(src[:,i]))
  if inv:
   for i,key in enumerate(['tailN0','tailN1','tailB0']):sim.write_across_shots(c.register_data[key],bits(want[:,i]))
  sim.do(c);errors=[]
  for i in range(N):
   if not np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(src[:,i])):errors.append(f'S{i}')
  if not inv:
   for i,key in enumerate(['tailN0','tailN1','tailB0']):
    if not np.array_equal(sim.read_across_shots(c.register_data[key]),bits(want[:,i])):errors.append(key)
  public={q for name,reg in c.register_data.items() if name!='endpoint_ghost' and (not inv or not name.startswith('tail')) for q in reg};junk=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in public]);r=dict(inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,errors=errors,scratch=bool(sim.read_across_shots(junk).any()),phase=bool(sim.read_phase_flipped_across_shots().any()));rows.append(r);print(json.dumps(r),flush=True);assert not errors and not r['scratch'] and not r['phase']
  if not inv:
   sim.do(ci)
   for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(src[:,i]))
   junk=km.array([km.q(q) for q in range(12*N,max(cf.num_qubits,ci.num_qubits))]);assert not sim.read_across_shots(junk).any();assert not sim.read_phase_flipped_across_shots().any()
 Path(f'results/raw_cf_endpoint_projective_n{N}_native.json').write_text(json.dumps(dict(status='PASS',shots=shots,rows=rows),indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);a=p.parse_args();test(a.n)
