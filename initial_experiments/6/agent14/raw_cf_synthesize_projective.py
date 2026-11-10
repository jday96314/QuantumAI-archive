"""Compose raw infinity endpoint, continued-fraction folds and normalization."""
import argparse,json,time
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import setup
from raw_cf_endpoint_projective import build as endpoint
from raw_cf_step import build as fold
from raw_cf_normalize import build as normalize
from schur_extract import body,mapped,BRE
from raw_cf_model import extract,synthesize
from quadratic_schur_model import F,node_maps


def layout(N):
 n=N//2;tail=[list(range(12*N+12*i,12*N+12*i+12)) for i in range(3)];A=tail[:2];B=tail[2:];free=[list(range(12*N+36+12*i,12*N+48+12*i)) for i in range(2)];steps=[]
 for i in reversed(range(n)):
  aa=list(range(24*i,24*i+12));bb=list(range(24*i+12,24*i+24));ports={'a':aa,'b':bb,**{f'N{j}':w for j,w in enumerate(A)},**{f'B{j}':w for j,w in enumerate(B)},**{f'Bt{j}':w for j,w in enumerate(free)}};steps.append(dict(i=i,size=len(A),ports=ports));A,B=B+free,A;free=[aa,bb]
 return dict(N=A,B=B,tail=tail,steps=steps,padding=free[0][:1],scratch=free[0][1:]+free[1]+list(range(12*N+60,12*N+240)))


def build(N,inverse=False):
 started=time.time();n=N//2;m,k,fw,bk,tr=setup();_,kap,_,_=node_maps(n);L=layout(N);b=km.CircuitBuilder()
 for i in range(N):b.create_quantum_register(12,f'S{i}')
 b.create_classical_register(49,'endpoint_ghost');b.create_classical_register(24*n,'head_ghost');path=Path(f'candidates/raw_cf_synthesize_projective_n{N}_{"inverse" if inverse else "forward"}.kmx');private=49+24*n
 with path.open('w') as f:
  f.write(str(b.finish_circuit())+'\n')
  def emit(c,ports,held=0,offset=0,scratch=None):
   qm={}
   for name,values in ports.items():
    assert len(c.register_data[name])==len(values),(name,len(c.register_data[name]),len(values));qm.update({q.id:v for q,v in zip(c.register_data[name],values)})
   missing=[i for i in range(c.num_qubits) if i not in qm];pool=list(range(12*N+60,12*N+240)) if scratch is None else scratch;available=[q for q in pool if q not in qm.values()];assert len(missing)<=len(available);qm.update(dict(zip(missing,available)))
   f.write(BRE.sub(lambda x:f'b{offset+int(x[1]) if int(x[1])<held else private+int(x[1])-held}',mapped(body(c),qm)))
  ep={f'S{i}':list(range(12*i,12*i+12)) for i in range(N)};ep.update({name:w for name,w in zip(['tailN0','tailN1','tailB0'],L['tail'])});ep_pool=list(range(12*N+36,12*N+240))
  norm={**{f'N{i}':w for i,w in enumerate(L['N'])},**{f'B{i}':w for i,w in enumerate(L['B'])},'padding':L['padding']}
  if inverse:emit(normalize(n,True),norm,scratch=L['scratch']+L['N'][-2]+L['N'][-1])
  else:emit(endpoint(N),ep,49,0,ep_pool)
  for st in (reversed(L['steps']) if inverse else L['steps']):emit(fold(st['size'],tr(kap[st['i']],fw),inverse),st['ports'],24,49+24*st['i'])
  if inverse:emit(endpoint(N,True),ep,49,0,ep_pool)
  else:emit(normalize(n),norm,scratch=L['scratch']+L['N'][-2]+L['N'][-1])
 c=km.Circuit.from_file(str(path),format='kmx');c.to_file(str(path.with_suffix('.kmb')),format='kmb');r=dict(N=N,inverse=inverse,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,seconds=time.time()-started);Path(f'results/raw_cf_synthesize_projective_n{N}_{"inverse" if inverse else "forward"}_build.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True);return c


def test(N=6,shots=256):
 m,_,fw,bk,tr=setup();n=N//2;rng=np.random.default_rng(142116+N);raw=[];wantA=[];wantB=[];norms=[];pads=[]
 for case in range(shots):
  d=case%(n+1);pos=rng.choice(np.arange(1,4096),d,replace=False).tolist();val=rng.integers(1,4096,d).tolist();heads,_=extract(F.syndrome(pos,val,N));A,B,exc,jet=synthesize(heads)
  if exc:continue
  pad=len(B)-1<n
  if pad:A=[0]+A;B=[0]+B
  lead=B[-1];inv=F.inv(lead);A=F.scale(A,inv);B=F.scale(B,inv);raw.append([tr(x,fw) for h in heads for x in h]);wantA.append([tr(x,fw) for x in A+[0]*(n+2-len(A))]);wantB.append([tr(x,fw) for x in B[:-1]]+[tr(lead,fw)]);pads.append(pad)
 raw=np.array(raw);wantA=np.array(wantA);wantB=np.array(wantB);shots=len(raw);cf=build(N);ci=build(N,True);L=layout(N);bits=lambda v,w=12:((np.asarray(v)[None,:]>>np.arange(w)[:,None])&1).astype(bool);sim=km.Simulator(batch_size=shots,seed=142116);sim.use_same_registers_as(cf)
 for i in range(N):sim.write_across_shots(cf.register_data[f'S{i}'],bits(raw[:,i]))
 sim.do(cf);errors=[]
 for key,bank,values in [('N',L['N'],wantA),('B',L['B'],wantB)]:
  for i,w in enumerate(bank):
   if not np.array_equal(sim.read_across_shots(km.array([km.q(q) for q in w])),bits(values[:,i])):errors.append(f'{key}{i}')
 assert np.array_equal(sim.read_across_shots(km.array([km.q(q) for q in L['padding']])),bits(pads,1));print(json.dumps(dict(stage='forward',errors=errors)),flush=True);assert not errors
 sim.do(ci)
 for i in range(N):assert np.array_equal(sim.read_across_shots(cf.register_data[f'S{i}']),bits(raw[:,i])),i
 scratch=km.array([km.q(q) for q in range(12*N,max(cf.num_qubits,ci.num_qubits))]);assert not sim.read_across_shots(scratch).any();assert not sim.read_across_shots(cf.register_data['endpoint_ghost']).any();assert not sim.read_across_shots(cf.register_data['head_ghost']).any();assert not sim.read_phase_flipped_across_shots().any();r=dict(status='PASS',N=N,shots=shots,forward=dict(Q=cf.num_qubits,T=cf.operation_statistics().num_toffoli_like),inverse=dict(Q=ci.num_qubits,T=ci.operation_statistics().num_toffoli_like));Path(f'results/raw_cf_synthesize_projective_n{N}_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);a=p.parse_args();test(a.n)
