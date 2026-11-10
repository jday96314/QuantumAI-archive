"""Completion-aware affine-chain synthesis, recycling the extraction frame.

The measured old frame has a pending private X outcome until the inverse
recreates it from original pivots. All extra fields/counters are explicit.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import kickmix as km
from extension_inverse import setup
from schur_affine_parameters_seedless import build as parameters
from schur_affine_dynamic_step_canonical import build as step
from schur_frame_product import build as frame_product
from schur_extract_repair import controls
from schur_extract_index import model,predicate
from schur_affine_outer_seedless import build as outer_step
from schur_extract import body,mapped,BRE
from schur_fold import reference
from quadratic_schur_model import F,node_maps


def layout(N):
 n=N//2;width=max(1,(n-1).bit_length());base=12*N+12+width;seed=list(range(12));free=[list(range(base,base+12)),list(range(base+12,base+24))];nodework=[base+24,base+25];scratch=list(range(12*N,12*N+12))+list(range(base+26,base+26+160));Nbank=[seed];Bbank=[];steps=[]
 for k in reversed(range(1,n)):
  ports=dict(u=list(range(24*k,24*k+12)),t=list(range(24*k+12,24*k+24)),Nlead=free[0],Btop=free[1]);ports.update({f'N{i}':w for i,w in enumerate(Nbank)});ports.update({f'B{i}':w for i,w in enumerate(Bbank)});steps.append(dict(k=k,d=len(Bbank),ports=ports));Nbank=Nbank+[free[0]];Bbank=Bbank+[free[1]];free=[ports['u'],ports['t']]
 ports=dict(b=list(range(12,24)),Nlead=free[0],Btop=free[1]);ports.update({f'N{i}':w for i,w in enumerate(Nbank)});ports.update({f'B{i}':w for i,w in enumerate(Bbank)});steps.append(dict(k=0,d=len(Bbank),ports=ports,outer=True))
 return dict(degree=[],seed=seed,steps=steps,N=Nbank,B=Bbank+[free[1]],zero_lead=free[0],index=list(range(12*N+12,12*N+12+width)),nodework=nodework,g=list(range(12*N,12*N+12)),scratch=scratch)


def classifier(n,k):
 b=km.CircuitBuilder();degree=b.create_quantum_register(n.bit_length(),'degree');c=b.create_quantum_register(1,'completed')[0];q=None
 b.flip_if_less_than(degree,k+1,target=c)
 x=b.finish_circuit();return km.Circuit(str(x)+'\n'+''.join(f'BIT_STORE0 b{i}\n' for i in range(x.num_bits)))


def build(N,reverse=False):
 begin=time.time();n=N//2;m,k,fw,bk,tr=setup();nodes,kap,_,_=node_maps(n);L=layout(N);cb=km.CircuitBuilder()
 for i in range(N):cb.create_quantum_register(12,f'S{i}')
 cb.create_quantum_register(12,'g');cb.create_quantum_register(len(L['index']),'index');cb.create_classical_register(12,'frame_ghost');cb.create_classical_register(n,'completion_ghost');header=str(cb.finish_circuit())+'\n';path=Path(f'candidates/schur_affine_synthesize_seedless_n{N}_{"inverse" if reverse else "forward"}.kmx')
 with path.open('w') as f:
  f.write(header)
  def emit(c,ports):
   qm={}
   for name,values in ports.items():
    reg=c.register_data[name];assert len(reg)==len(values),(name,len(reg),len(values));qm.update({q.id:v for q,v in zip(reg,values)})
   missing=[i for i in range(c.num_qubits) if i not in qm];available=[q for q in L['scratch'] if q not in qm.values()];assert len(missing)<=len(available);qm.update(dict(zip(missing,available)));offset=12 if 'completion_ghost' in c.register_data else 12+n;f.write(BRE.sub(lambda x:f'b{offset+int(x[1])}',mapped(body(c),qm)))
  pp={f'S{i}':list(range(12*i,12*i+12)) for i in range(N)};
  if not reverse:
   for bit,q in enumerate(L['g']):f.write(f'BIT_STORE0 b{bit}\nHMR q{q} b{bit}\n')
   emit(parameters(N),pp)
  for s in (reversed(L['steps']) if reverse else L['steps']):
   k=s['k'];cc=controls(N,k);deltas=[tr(F.mul(de,de)^de,fw) for _,de in cc];work=L['nodework'][:len(cc)];ports={**s['ports'],'flags':work}
   for (flag,_),target in zip(cc,work):emit(predicate(len(L['index']),flag+1),dict(index=L['index'],target=[target]))
   emit((outer_step if s.get('outer') else step)(s['d'],tr(kap[k],fw),deltas,reverse),ports)
   for (flag,_),target in reversed(list(zip(cc,work))):emit(predicate(len(L['index']),flag+1),dict(index=L['index'],target=[target]))
  if reverse:
   emit(parameters(N,True),pp);gp={f'S{i}':list(range(12*i,12*i+12)) for i in range(N)};gp['g']=L['g'];emit(frame_product(N),gp)
   for bit,q in enumerate(L['g']):f.write(f'Z q{q} if b{bit}\nBIT_STORE0 b{bit}\n')
 c=km.Circuit.from_file(str(path),format='kmx');c.to_file(str(path.with_suffix('.kmb')),format='kmb');r=dict(N=N,reverse=reverse,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,seconds=time.time()-begin,numerator_qubits=L['N'],denominator_qubits=L['B'],degree_qubits=L['degree'],scope='complete affine synthesis pair, pending measured-frame phase forward; native/full verification required');Path(f'results/schur_affine_synthesize_seedless_n{N}_{"inverse" if reverse else "forward"}_build.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if not isinstance(v,list)}),flush=True);return c


def test(N=6,shots=256):
 m,k,fw,bk,tr=setup();n=N//2;rng=np.random.default_rng(142067+N);_,kap,_,_=node_maps(n);source=[];flags=[];wn=[];wb=[];gg=[];ff=[];dd=[];actual=0
 for case in range(shots):
  if case%2==0:
   d=(case//2)%(n+1);pos=rng.choice(np.arange(1,4096),d,replace=False).tolist();val=rng.integers(1,4096,d).tolist();vv,g,f,_=model(F.syndrome(pos,val,N));v=[tr(x,fw) for x in vv]
   if any(not v[2*i+1] for i in range(d)):continue
   actual+=1
  else:
   d=case%(n+1);v=rng.integers(0,4096,N).tolist();v[1::2]=rng.integers(1,4096,n).tolist();v[2*d:]=[0]*(N-2*d);f=[]
   for i in range(n-1):f.append(int(i<d-1 and not any(f) and rng.random()<.3))
  nn=[];bb=[];g=1
  for i in range(n):g=m.gmul(m.gmul(v[2*i+1] or 1,v[2*i+1] or 1),m.ginv(g))
  for i in reversed(range(n)):
   kk=kap[i]
   for fl,de in controls(N,i):
    if f[fl]:kk^=F.mul(de,de)^de
   nn,bb=reference(v[2*i],v[2*i+1],nn,bb,tr(kk,fw))
  source.append(v);flags.append(next((i+1 for i,z in enumerate(f) if z),0));wn.append(nn);wb.append(bb);gg.append(g);ff.append(m.ginv(v[2*d-1]) if d else 1);dd.append(d)
 src=np.array(source);flags=((np.array(flags)[None,:]>>np.arange(len(layout(N)['index']))[:,None])&1).astype(bool);wn=np.array(wn);wb=np.array(wb);gg=np.array(gg);ff=np.array(ff);dd=np.array(dd);shots=len(src);cf=build(N);ci=build(N,True);L=layout(N);bits=lambda v,w=12:((v[None,:]>>np.arange(w)[:,None])&1).astype(bool);reports=[]
 for inv,c in [(False,cf),(True,ci)]:
  sim=km.Simulator(batch_size=shots,seed=142067+inv);sim.use_same_registers_as(c);sim.write_across_shots(c.register_data['g'],bits(np.zeros(shots,dtype=int) if inv else gg));sim.write_across_shots(c.register_data['index'],flags)
  if inv:
   for bank,values in [(L['N'],wn),(L['B'],wb)]:
    for i,w in enumerate(bank):sim.write_across_shots(km.array([km.q(q) for q in w]),bits(values[:,i]))

  else:
   for i in range(N):sim.write_across_shots(c.register_data[f'S{i}'],bits(src[:,i]))
  sim.do(c);errors={}
  if inv:
   for i in range(N):
    bad=np.flatnonzero((sim.read_across_shots(c.register_data[f'S{i}'])!=bits(src[:,i])).any(0))
    if len(bad):errors[f'S{i}']=bad[:8].tolist()
   used=set(range(12*N+12+len(L['index'])));gout=gg
  else:
   for name,bank,values in [('N',L['N'],wn),('B',L['B'],wb)]:
    for i,w in enumerate(bank):
     bad=np.flatnonzero((sim.read_across_shots(km.array([km.q(q) for q in w]))!=bits(values[:,i])).any(0))
     if len(bad):errors[f'{name}{i}']=bad[:8].tolist()
   used={q for w in L['N']+L['B'] for q in w}|set(L['degree'])|set(L['g'])|set(L['index']);gout=np.zeros(shots,dtype=int)
  bad=np.flatnonzero((sim.read_across_shots(c.register_data['g'])!=bits(gout)).any(0))
  if len(bad):errors['g']=bad[:8].tolist()
  assert np.array_equal(sim.read_across_shots(c.register_data['index']),flags);scratch=km.array([km.q(i) for i in range(c.num_qubits) if i not in used]);junk=np.flatnonzero(sim.read_across_shots(scratch).any(0));mask=sim.read_across_shots(c.register_data['frame_ghost']);expect=np.zeros(shots,dtype=bool) if inv else ((mask&bits(gg)).sum(0)%2).astype(bool);cmask=sim.read_across_shots(c.register_data['completion_ghost']);expect^=((cmask&(np.arange(n)[:,None]>=dd[None,:])).sum(0)%2).astype(bool) if not inv else False;phase=np.flatnonzero(sim.read_phase_flipped_across_shots()!=expect);r=dict(inverse=inv,errors=errors,junk=junk[:8].tolist(),phase=phase[:8].tolist());reports.append(r);print(json.dumps(r),flush=True);assert not errors and not len(junk) and not len(phase)
  if not inv:
   sim.do(ci)
   for i in range(N):assert np.array_equal(sim.read_across_shots(c.register_data[f'S{i}']),bits(src[:,i])),i
   assert np.array_equal(sim.read_across_shots(c.register_data['g']),bits(gg));assert not sim.read_across_shots(c.register_data['frame_ghost']).any();assert not sim.read_across_shots(c.register_data['completion_ghost']).any();assert not sim.read_phase_flipped_across_shots().any();pair_scratch=km.array([km.q(q) for q in range(max(cf.num_qubits,ci.num_qubits)) if q>=12*N+12+len(L['index'])]);assert not sim.read_across_shots(pair_scratch).any()
 result=dict(status='PASS',N=N,shots=shots,true_syndrome_cases=actual,rows=reports,scope='all prefix lengths, generic pivots and actual repaired syndromes, variable node descriptors, exact coefficient model, degree/boundary/frame, forward ghost phase, fresh inverse and paired cleanup');Path(f'results/schur_affine_synthesize_seedless_n{N}_native.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=6);p.add_argument('--shots',type=int,default=256);a=p.parse_args();test(a.n,a.shots)
