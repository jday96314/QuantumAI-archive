"""A source-resident three-field endpoint with monic final cofactor gauge."""
import argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import Adapter,setup
from fast_dqi import FastBuilder
from lean_emitter import create
from cubic_low_workspace_inverse_compact import bind
from raw_three_jet_step import build as jet_step
from raw_three_jet_model import endpoint
from raw_cf_model import extract
from quadratic_schur_model import F
from schur_extract import body,mapped


def boundary(inverse=False):
 m=setup()[0];m.inv_word=bind(m);b=FastBuilder(km.CircuitBuilder());D=b.create_quantum_register(12,'d');P=b.create_quantum_register(12,'p');Q=b.create_quantum_register(12,'q');ab=Adapter(b);sc,E=create(m,ab);d=[q.id for q in D];q=[v.id for v in Q];z=b.alloc_qubits(1)[0]
 def detect():
  with m.eq_test(E,q,0) as zero:E.cx(zero,z.id)
 detect()
 if inverse:b.cswap(z,D,P)
 b.cx(z,Q[0])
 with m.inv_word(E,q) as qi,m.opened(E,'three_jet_boundary_product',12) as tmp:m.mon_slot(E,d,q if inverse else qi,qi if inverse else q,tmp)
 b.cx(z,Q[0])
 if not inverse:b.cswap(z,D,P)
 detect();b.free(km.array([z]));ab.release_recycle('three_boundary_emitter',sc);c=b.finish_circuit();return km.Circuit(str(c)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(c.num_bits)))


def build(N,inverse=False):
 n=N//2;b=km.CircuitBuilder()
 for i in range(N):b.create_quantum_register(12,f'S{i}')
 tail=[list(range(12)),list(range(12*N,12*N+12)),list(range(12,24))];headers=str(b.finish_circuit())+'\n'+''.join(''.join(f'APPEND_TO_REGISTER q{q} r{N+i}\n' for q in w)+f'REGISTER r{N+i} "{name}"\n' for i,(name,w) in enumerate(zip(['tailN0','tailN1','tailB0'],tail)));out=[headers];pool=list(range(12*N+12,12*N+160))
 def gate(op,*qs):out.append(op+' '+ ' '.join(f'q{q}' for q in qs)+'\n')
 def emit(c,ports):
  qm={q.id:v for name,w in ports.items() for q,v in zip(c.register_data[name],w)};qm.update(dict(zip([q for q in range(c.num_qubits) if q not in qm],pool)));out.append(mapped(body(c),qm))
 base=dict(d=tail[0],p=tail[1],q=tail[2])
 if inverse:emit(boundary(True),base)
 else:
  for q,r in zip(tail[2],tail[0]):gate('CX',q,r)
  gate('X',tail[1][0])
 c=jet_step(inverse)
 for i in (reversed(range(1,n)) if inverse else range(1,n)):emit(c,dict(a=list(range(24*i,24*i+12)),b=list(range(24*i+12,24*i+24)),**base))
 if inverse:
  gate('X',tail[1][0])
  for q,r in zip(tail[2],tail[0]):gate('CX',q,r)
 else:emit(boundary(),base)
 return km.Circuit(''.join(out))


def test(N=6,shots=512):
 m,_,fw,bk,tr=setup();n=N//2;rng=np.random.default_rng(142138+N);ss=[];want=[]
 for case in range(shots):
  d=case%(n+1);pos=rng.choice(np.arange(1,4096),d,replace=False).tolist();val=rng.integers(1,4096,d).tolist();hh,_=extract(F.syndrome(pos,val,N));A,B,e,j=endpoint(hh)
  if e:continue
  A=A+[0]*(2-len(A));src=[tr(v,fw) for h in hh for v in h];dst=src.copy();dst[0]=tr(A[0],fw);dst[1]=tr(B[0] if B else 0,fw);ss.append(src);want.append(dst+[tr(A[1],fw)])
 ss=np.array(ss);want=np.array(want);shots=len(ss);bits=lambda v:((np.asarray(v)[None,:]>>np.arange(12)[:,None])&1).astype(bool);cf=build(N);ci=build(N,True);rows=[]
 for inv,c in [(False,cf),(True,ci)]:
  sim=km.Simulator(batch_size=shots,seed=142138+inv);sim.use_same_registers_as(c)
  for i in range(N):sim.write_across_shots(c.register_data[f'S{i}'],bits(want[:,i] if inv else ss[:,i]))
  if inv:sim.write_across_shots(c.register_data['tailN1'],bits(want[:,-1]))
  sim.do(c)
  for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(ss[:,i] if inv else want[:,i])),i
  assert np.array_equal(sim.read_across_shots(c.register_data['tailN1']),bits(np.zeros(shots,dtype=int) if inv else want[:,-1]));public={q for name,reg in c.register_data.items() if name.startswith('S') or (not inv and name=='tailN1') for q in reg};junk=km.array([km.q(q) for q in range(c.num_qubits) if km.q(q) not in public]);assert not sim.read_across_shots(junk).any();assert not sim.read_phase_flipped_across_shots().any();rows.append(dict(inverse=inv,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like));print(json.dumps(rows[-1]),flush=True)
 r=dict(status='PASS',N=N,shots=shots,rows=rows);Path(f'results/raw_three_jet_endpoint_n{N}_native.json').write_text(json.dumps(r,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);a=p.parse_args();test(a.n)
