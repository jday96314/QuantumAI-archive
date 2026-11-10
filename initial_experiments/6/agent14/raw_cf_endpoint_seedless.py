"""Use the first source head as the normalized terminal jet workspace."""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from projective_head_jet import build as first_step
from projective_jet_step import build as jet_step
from raw_cf_model import extract,endpoint
from schur_extract import body,mapped,BRE
from extension_inverse import setup
from quadratic_schur_model import F


def build(N,inverse=False):
 n=N//2;b=km.CircuitBuilder()
 for i in range(N):b.create_quantum_register(12,f'S{i}')
 tail=[list(range(12)),list(range(12,24)),list(range(12*N,12*N+12))];tag=12*N+12
 headers=str(b.finish_circuit())+'\n'+''.join(''.join(f'APPEND_TO_REGISTER q{q} r{N+i}\n' for q in w)+f'REGISTER r{N+i} "{name}"\n' for i,(name,w) in enumerate(zip(['tailN0','tailN1','tailB0'],tail)));out=[headers];pool=tail[2]+list(range(tag+1,tag+90))
 def gate(op,*qs):out.append(op+' '+ ' '.join(f'q{q}' for q in qs)+'\n')
 def emit(c,ports):
  qm={q.id:v for name,w in ports.items() for q,v in zip(c.register_data[name],w)};qm.update(dict(zip([q for q in range(c.num_qubits) if q not in qm],pool)));out.append(mapped(body(c),qm))
 def boundary():
  if inverse:gate('X',tag);gate('CX',tail[2][0],tag);gate('X',tail[2][0]);gate('CX',tag,tail[2][0]);gate('CX',tag,tail[0][0])
  else:gate('CX',tag,tail[0][0]);gate('X',tail[2][0]);gate('CX',tag,tail[2][0]);gate('CX',tail[2][0],tag);gate('X',tag)
 base={'h':tail[1],'k':tail[0],'infinity':[tag]}
 if inverse:boundary()
 else:emit(first_step(),base)
 c=jet_step(inverse)
 for i in (reversed(range(1,n)) if inverse else range(1,n)):emit(c,dict(a=list(range(24*i,24*i+12)),b=list(range(24*i+12,24*i+24)),**base))
 if inverse:emit(first_step(True),base)
 else:boundary()
 return km.Circuit(''.join(out))


def test(N=6,shots=512):
 m,_,fw,bk,tr=setup();n=N//2;rng=np.random.default_rng(142123+N);ss=[];want=[]
 for case in range(shots):
  d=case%(n+1);pos=rng.choice(np.arange(1,4096),d,replace=False).tolist();val=rng.integers(1,4096,d).tolist();hh,_=extract(F.syndrome(pos,val,N));a,b,e,j=endpoint(hh)
  if e:continue
  a=a+[0]*(2-len(a));src=[tr(v,fw) for h in hh for v in h];dst=src.copy();dst[:2]=[tr(a[0],fw),tr(a[1],fw)];ss.append(src);want.append(dst+[int(bool(b))])
 ss=np.array(ss);want=np.array(want);shots=len(ss);bits=lambda v:((np.asarray(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool);cf=build(N);ci=build(N,True);rows=[]
 for inv,c in [(False,cf),(True,ci)]:
  sim=km.Simulator(batch_size=shots,seed=142123+inv);sim.use_same_registers_as(c)
  for i in range(N):sim.write_across_shots(c.register_data[f'S{i}'],bits(want[:,i] if inv else ss[:,i]))
  if inv:sim.write_across_shots(c.register_data['tailB0'],bits(want[:,-1]))
  sim.do(c)
  for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(ss[:,i] if inv else want[:,i])),i
  assert np.array_equal(sim.read_across_shots(c.register_data['tailB0']),bits(np.zeros(shots,dtype=int) if inv else want[:,-1]))
  public={q for name,reg in c.register_data.items() if name.startswith('S') or (not inv and name=='tailB0') for q in reg};junk=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in public]);assert not sim.read_across_shots(junk).any();assert not sim.read_phase_flipped_across_shots().any();r=dict(inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like);rows.append(r);print(json.dumps(r),flush=True)
 Path(f'results/raw_cf_endpoint_seedless_n{N}_native.json').write_text(json.dumps(dict(status='PASS',N=N,shots=shots,rows=rows),indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);a=p.parse_args();test(a.n)
