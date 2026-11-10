"""Narrow prefix plus leading-gap suffix, with shared quotient synthesis."""
from pathlib import Path
import argparse,json
import prefix_counter_oracle as base
import prefix_high_copy_jacobi as prefix
import lookahead_jacobi as look
import compact_lookahead_parameters as compact
import mixed_quotient_jacobi as mixed
from import_published import load
from projective_polynomial_inverse import bind
from compact_blocked_cfraction import Arithmetic
import dqi_budget_probe as probe
HERE=Path(__file__).resolve().parent

def extract(m,E,math,W,counts,through,inverse=False):
 stages=[];k=0
 for b in counts:stages.append((k,b));k+=b
 for k,b in (reversed(stages) if inverse else stages):
  if b==1:prefix.pair(m,E,math,W[2*k:],inverse)
  elif k<through:prefix.block(m,E,math,W[2*k:],b,inverse)
  else:
   assert b==2
   look.block(m,E,math,W[2*k:],inverse)

def emit(N,pairs,through,suffix,inverse_limit,profile=False):
 assert through%2==0 and suffix==2
 look.parameters=compact.parameters;m=load();m.inv_word=bind(m,limit=inverse_limit)
 if profile:
  Base=m.Em
  class Counter(Base):
   def emit_mac(self,rt):self.b._cur_stats['toffoli']+=sum(t[0]=='P' for t in rt)
   def emit_free(self,rt):pass
  m.Em=Counter
 b=m.Builder()
 if profile:b._items=probe.count_module().Builder()._items
 W=[b.qubits(12) for _ in range(N)];E=m.Em(b,*b.qubits(3),m.Alloc(b));math=Arithmetic(m,E);stages=[];counts=prefix.schedule(N//2,pairs,through,suffix)
 def stage(name,fn):
  before=b._cur_stats['toffoli'];ret=fn();stages.append(dict(name=name,T=b._cur_stats['toffoli']-before,peak=b.peak_qubits));return ret
 stage('extract',lambda:extract(m,E,math,W,counts,through))
 n=N//2;ng=sum(k>=through for k in range(0,n-1,2))
 with m.opened(E,'mixed_terminal',12) as c,m.opened(E,'mixed_metadata',n.bit_length()+1) as metadata,m.opened(E,'mixed_gap_flags',ng) as gaps:
  F,G=stage('synthesize',lambda:mixed.synthesize(m,E,W,c,metadata,gaps,through=through));hook=len(b._items);held=[*c,*metadata,*gaps]
  stage('synthesize_inverse',lambda:mixed.synthesize(m,E,W,c,metadata,gaps,True,through))
 stage('extract_inverse',lambda:extract(m,E,math,W,counts,through,True))
 return m,b,W,F,G,held,hook,dict(N=N,prefix=pairs,through=through,suffix=suffix,inverse_limit=inverse_limit,counts=counts,core_Q=b.peak_qubits,core_T=b._cur_stats['toffoli'],stages=stages,scope='narrow regular prefix and leading-gap suffix, shared degree1/2 quotient reconstruction; full verification required')

def build(N=14,M=127,pairs=2,through=2,suffix=2,inverse_limit=935,tag=None):
 base.emit=emit
 return base.build(N,M,pairs,through,suffix,inverse_limit,tag=tag or f'mixed_prefix_n{N}_m{M}_b{pairs}_t{through}_l{inverse_limit}')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--profile',action='store_true');p.add_argument('--n',type=int,default=70);p.add_argument('--m',type=int,default=4095);p.add_argument('--pairs',type=int,default=6);p.add_argument('--through',type=int,default=18);p.add_argument('--inverse-limit',type=int,default=935);a=p.parse_args()
 if a.profile:
  *_,r=emit(a.n,a.pairs,a.through,2,a.inverse_limit,True);r.update(projected_scan_T=80*a.m,projected_total_T=r['core_T']+80*a.m,status='PROFILE_ONLY_NO_FULL_NATIVE_ACCEPTANCE');(HERE/'results'/f'mixed_prefix_b{a.pairs}_t{a.through}_l{a.inverse_limit}_profile.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
 else:build(a.n,a.m,a.pairs,a.through,2,a.inverse_limit)
