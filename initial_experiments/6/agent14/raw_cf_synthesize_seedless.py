"""Source-resident terminal jet and proper outer contraction."""
import argparse,json,time
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from cubic_low_workspace_inverse import bind
from fast_dqi import FastBuilder
from lean_emitter import create
from raw_cf_endpoint_seedless import build as endpoint
from raw_cf_step import build as fold
from raw_cf_outer_seedless import build as outer
from schur_extract import body,mapped,BRE
from raw_cf_model import extract,synthesize
from quadratic_schur_model import F,node_maps


def layout(N):
 n=N//2;tail=[list(range(12)),list(range(12,24)),list(range(12*N,12*N+12))];A=tail[:2];B=tail[2:];free=[list(range(12*N+12+12*i,12*N+24+12*i)) for i in range(2)];steps=[]
 for i in reversed(range(1,n)):
  aa=list(range(24*i,24*i+12));bb=list(range(24*i+12,24*i+24));ports={'a':aa,'b':bb,**{f'N{j}':w for j,w in enumerate(A)},**{f'B{j}':w for j,w in enumerate(B)},**{f'Bt{j}':w for j,w in enumerate(free)}};steps.append(dict(i=i,size=len(A),ports=ports));A,B=B+free,A;free=[aa,bb]
 pad=[12*N+36];return dict(N=B,B=A,tail=tail,steps=steps,preN=A,preB=B,padding=pad,free=free,scratch=free[0]+free[1]+list(range(12*N+37,12*N+200)))


def normalize(n,inverse=False):
 m=setup()[0];m.inv_word=bind(m);b=FastBuilder(km.CircuitBuilder());N=[b.create_quantum_register(12,f'N{i}') for i in range(n)];B=[b.create_quantum_register(12,f'B{i}') for i in range(n+1)];ab=Adapter(b);sc,E=create(m,ab);lead=[q.id for q in B[-1]]
 with m.inv_word(E,lead) as li,m.opened(E,'seedless_normalize_product',12) as tmp:
  for w in N+B[:n]:m.mon_slot(E,[q.id for q in w],lead if inverse else li,li if inverse else lead,tmp)
 ab.release_recycle('seedless_normalize_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def build(N,inverse=False):
 started=time.time();n=N//2;m,k,fw,bk,tr=setup();_,kap,_,_=node_maps(n);L=layout(N);b=km.CircuitBuilder()
 for i in range(N):b.create_quantum_register(12,f'S{i}')
 b.create_classical_register(24*n,'head_ghost');path=Path(f'candidates/raw_cf_synthesize_seedless_n{N}_{"inverse" if inverse else "forward"}.kmx');private=24*n
 with path.open('w') as f:
  f.write(str(b.finish_circuit())+'\n')
  def emit(c,ports,held=0,offset=0,scratch=None):
   qm={q.id:v for name,values in ports.items() for q,v in zip(c.register_data[name],values)};missing=[i for i in range(c.num_qubits) if i not in qm];pool=list(range(12*N+36,12*N+200)) if scratch is None else scratch;available=[q for q in pool if q not in qm.values()];assert len(missing)<=len(available);qm.update(dict(zip(missing,available)));f.write(BRE.sub(lambda x:f'b{offset+int(x[1]) if int(x[1])<held else private+int(x[1])-held}',mapped(body(c),qm)))
  ep={f'S{i}':list(range(12*i,12*i+12)) for i in range(N)};ep.update({name:w for name,w in zip(['tailN0','tailN1','tailB0'],L['tail'])});ep_pool=list(range(12*N+12,12*N+200));norm={**{f'N{i}':w for i,w in enumerate(L['N'])},**{f'B{i}':w for i,w in enumerate(L['B'])}};op={**{f'N{i}':w for i,w in enumerate(L['preN'])},**{f'B{i}':w for i,w in enumerate(L['preB'])}}
  def padding(back=False):
   flag=L['padding'][0]
   def detect():
    reg=L['preN'][-1];cb=FastBuilder(km.CircuitBuilder());rr=cb.create_quantum_register(12,'lead');fl=cb.create_quantum_register(1,'padding')[0];cb.flip_if_equal(rr,0,target=fl);emit(cb.finish_circuit(),{'lead':reg,'padding':[flag]},scratch=L['scratch'])
   if not back:detect()
   for bank in [L['preN'],L['preB']]:
    for i in (range(len(bank)-1) if back else reversed(range(len(bank)-1))):
     for q,r in zip(bank[i],bank[i+1]):f.write(f'CSWAP q{flag} q{q} q{r}\n')
   if back:detect()
  if inverse:
   emit(normalize(n,True),norm,scratch=L['scratch']);emit(outer(n,tr(kap[0],fw),True),op,24,0,scratch=L['scratch']);padding(True)
  else:emit(endpoint(N),ep,scratch=ep_pool)
  for st in (reversed(L['steps']) if inverse else L['steps']):emit(fold(st['size'],tr(kap[st['i']],fw),inverse),st['ports'],24,24*st['i'])
  if inverse:emit(endpoint(N,True),ep,scratch=ep_pool)
  else:
   padding();emit(outer(n,tr(kap[0],fw)),op,24,0,scratch=L['scratch']);emit(normalize(n),norm,scratch=L['scratch'])
 c=km.Circuit.from_file(str(path),format='kmx');c.to_file(str(path.with_suffix('.kmb')),format='kmb');r=dict(N=N,inverse=inverse,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,seconds=time.time()-started);Path(f'results/raw_cf_synthesize_seedless_n{N}_{"inverse" if inverse else "forward"}_build.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True);return c


def test(N=6,shots=256):
 m,_,fw,bk,tr=setup();n=N//2;rng=np.random.default_rng(142125+N);raw=[];wantA=[];wantB=[];pads=[]
 for case in range(shots):
  d=case%(n+1);pos=rng.choice(np.arange(1,4096),d,replace=False).tolist();val=rng.integers(1,4096,d).tolist();heads,_=extract(F.syndrome(pos,val,N));A,B,exc,jet=synthesize(heads)
  if exc:continue
  pad=len(B)-1<n
  if pad:A=[0]+A;B=[0]+B
  lead=B[-1];inv=F.inv(lead);A=F.scale(A,inv);B=F.scale(B,inv);raw.append([tr(x,fw) for h in heads for x in h]);wantA.append([tr(x,fw) for x in A+[0]*(n-len(A))]);wantB.append([tr(x,fw) for x in B[:-1]]+[tr(lead,fw)]);pads.append(pad)
 raw=np.array(raw);wantA=np.array(wantA);wantB=np.array(wantB);shots=len(raw);cf=build(N);ci=build(N,True);L=layout(N);bits=lambda v,w=12:((np.asarray(v)[None,:]>>np.arange(w)[:,None])&1).astype(bool);sim=km.Simulator(batch_size=shots,seed=142125);sim.use_same_registers_as(ci)
 for i in range(N):sim.write_across_shots(cf.register_data[f'S{i}'],bits(raw[:,i]))
 sim.do(cf);errors=[]
 for key,bank,values in [('N',L['N'],wantA),('B',L['B'],wantB)]:
  for i,w in enumerate(bank):
   if not np.array_equal(sim.read_across_shots(km.array([km.q(q) for q in w])),bits(values[:,i])):errors.append(f'{key}{i}')
 assert np.array_equal(sim.read_across_shots(km.array([km.q(q) for q in L['padding']])),bits(pads,1));print(json.dumps(dict(stage='forward',errors=errors)),flush=True);assert not errors
 sim.do(ci)
 for i in range(N):assert np.array_equal(sim.read_across_shots(cf.register_data[f'S{i}']),bits(raw[:,i])),i
 scratch=km.array([km.q(q) for q in range(12*N,max(cf.num_qubits,ci.num_qubits))]);assert not sim.read_across_shots(scratch).any();assert not sim.read_across_shots(cf.register_data['head_ghost']).any();assert not sim.read_phase_flipped_across_shots().any();r=dict(status='PASS',N=N,shots=shots,forward=dict(Q=cf.num_qubits,T=cf.operation_statistics().num_toffoli_like),inverse=dict(Q=ci.num_qubits,T=ci.operation_statistics().num_toffoli_like));Path(f'results/raw_cf_synthesize_seedless_n{N}_native.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);a=p.parse_args();test(a.n)
