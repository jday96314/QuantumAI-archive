"""Complete research oracle: prefix Jacobi blocks and compact completion.

The point body is the Apache-2.0 round-3 recursive Forney implementation.
All prior artifacts and host verifiers are read-only. This builder does not
assert decoding correctness on singular prefix or interior-gap charts.
"""
from pathlib import Path
import argparse,hashlib,json,subprocess,time
import numpy as np
import kickmix as km
from import_published import load
from projective_polynomial_inverse import bind
from compact_blocked_cfraction import Arithmetic
import prefix_high_copy_jacobi as prefix
import counter_quotient_jacobi as quotient
from published_dqi import isomorphism
from fast_dqi import FastBuilder
from recursive_forney import RecursiveForney
from tower_inverse import map_ops,apply_map,linear_xor

HERE=Path(__file__).resolve().parent

def emit(N,pairs,through,suffix,inverse_limit):
 m=load();m.inv_word=bind(m,limit=inverse_limit);b=m.Builder();W=[b.qubits(12) for _ in range(N)]
 E=m.Em(b,*b.qubits(3),m.Alloc(b));math=Arithmetic(m,E);stages=[]
 counts=prefix.schedule(N//2,pairs,through,suffix)
 def stage(name,fn):
  before=b._cur_stats['toffoli'];ret=fn();stages.append(dict(name=name,T=b._cur_stats['toffoli']-before,peak=b.peak_qubits));return ret
 stage('extract',lambda:prefix.extract(m,E,math,W,counts=counts,through=through))
 with m.opened(E,'quotient_terminal',12) as c,m.opened(E,'quotient_metadata',(N//2).bit_length()+1) as metadata:
  F,G=stage('synthesize',lambda:quotient.synthesize(m,E,W,c,metadata));hook=len(b._items)
  held=[*c,*metadata]
  stage('synthesize_inverse',lambda:quotient.synthesize(m,E,W,c,metadata,True))
 stage('extract_inverse',lambda:prefix.extract(m,E,math,W,True,counts,through))
 return m,b,W,F,G,held,hook,dict(N=N,prefix=pairs,through=through,suffix=suffix,inverse_limit=inverse_limit,counts=counts,core_Q=b.peak_qubits,core_T=b._cur_stats['toffoli'],stages=stages)

def build(N=14,M=127,pairs=3,through=3,suffix=2,inverse_limit=935,tag=None):
 started=time.time();m,b,W,F,G,held,hook,report=emit(N,pairs,through,suffix,inverse_limit)
 tag=tag or f'prefix_counter_n{N}_m{M}_b{pairs}_t{through}_s{suffix}_l{inverse_limit}'
 stem=HERE/'candidates'/tag;tape=HERE/'parameters'/f'{tag}.npy'
 rows=np.asarray(b._items._rows,dtype=np.uint32).reshape(-1,9);np.save(tape,rows)
 kb=FastBuilder(km.CircuitBuilder());S=[kb.create_quantum_register(12,f'S{i}') for i in range(N)];masks=[kb.create_classical_register(12,f'E{i}') for i in range(M)]
 q=[x for w in S for x in w]+list(kb.alloc_qubits(b.next_qubit-12*N));fw,bk,tr=isomorphism(m)
 for word in S:apply_map(kb,word,map_ops(fw))
 kb.info_tag('PREFIX_SCAN_BEGIN');keep=set(held);kb.free(km.array([q[i] for i in range(12*N,b.next_qubit) if i not in keep]))
 F=[km.array([q[i] for i in word]) for word in F];G=[km.array([q[i] for i in word]) for word in G]
 for word in F+G:apply_map(kb,word,map_ops(bk))
 gf=km.GF2Field(12)
 class Forney(RecursiveForney):gf=km.GF2Field(12)
 body=Forney();derivative,root,numerator=F[1],F[0],G[0]
 def lin(word,target,c):
  if c:linear_xor(kb,word,target,[gf.mul(c,1<<i) for i in range(12)])
 for j,mask in enumerate(masks):
  a=gf.pow(2,j%4095)
  def d():
   for i in range(3,len(F),2):lin(F[i],derivative,gf.pow(a,i-1))
  def r():
   lin(derivative,root,a)
   for i in range(2,len(F),2):lin(F[i],root,gf.pow(a,i))
  def p():
   for i in range(1,len(G)):lin(G[i],numerator,gf.pow(a,i))
  # e=G(a)/(a F'(a)); apply the dual of multiplication by 1/a
  # to the run-time classical mask and restore it immediately afterward.
  scale=gf.pow(a,4094);maskops=map_ops([sum(((gf.mul(scale,1<<z)>>i)&1)<<z for z in range(12)) for i in range(12)])
  d();r();p();apply_map(kb,mask,maskops);body._forney_kernel(kb,root,derivative,numerator,mask);apply_map(kb,mask,maskops,True);p();r();d()
  if M>127 and j%512==0:print(json.dumps(dict(stage='scan',j=j,seconds=time.time()-started)),flush=True)
 for word in F+G:apply_map(kb,word,map_ops(bk),True)
 # Some native measurement helpers free their classical outcome without
 # zeroing it. Clear the private scan outcomes before another invocation.
 # No public mask or quantum input is measured or changed by this cleanup.
 scan_bits=kb.finish_circuit().num_bits
 kb.bit_store(km.array([km.b(i) for i in range(12*M,scan_bits)]),False)
 kb.info_tag('PREFIX_SCAN_END')
 for word in S:apply_map(kb,word,map_ops(fw),True)
 wrapper=stem.with_suffix('.wrapper.kmx');kb.finish_circuit().to_file(str(wrapper),format='kmx')
 parts=[stem.with_suffix('.'+part+'.kmx') for part in ('pre','scan','post')];handles=[p.open('w') for p in parts];part=0
 for line in wrapper.read_text().splitlines(True):
  if 'PREFIX_SCAN_BEGIN' in line:part=1
  elif 'PREFIX_SCAN_END' in line:part=2
  else:handles[part].write(line)
 for h in handles:h.close()
 empty=HERE/'parameters/empty';empty.write_text('');out=stem.with_suffix('.kmx')
 subprocess.run([str(HERE/'../../4/codex/translate_prefix'),str(tape),str(hook),str(12*M),*[str(p) for p in parts],str(out),str(b.next_out),str(empty),str(len(rows))],check=True,capture_output=True)
 circuit=km.Circuit.from_file(str(out),format='kmx');path=stem.with_suffix('.kmb');circuit.to_file(str(path),format='kmb')
 report.update(M=M,Q=circuit.num_qubits,T=circuit.operation_statistics().num_toffoli_like,path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),seconds=time.time()-started,status='BUILT_REQUIRES_NATIVE_VERIFICATION',scope='complete approximate-chart research oracle; no submission acceptance')
 assert report['T']==report['core_T']+80*M
 (HERE/'results'/f'{tag}_build.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
 return circuit,report

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=14);p.add_argument('--m',type=int,default=127);p.add_argument('--pairs',type=int,default=3);p.add_argument('--through',type=int,default=3);p.add_argument('--suffix',type=int,default=2);p.add_argument('--inverse-limit',type=int,default=935);a=p.parse_args();build(a.n,a.m,a.pairs,a.through,a.suffix,a.inverse_limit)
