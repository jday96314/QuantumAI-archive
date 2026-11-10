"""Full arithmetic profiles for narrow prefix extraction and count completion.

This emits every arithmetic stage into the verified count-only emitter. The
80T point body is added explicitly, separately from the emitted arithmetic.
No native full-task acceptance or decoding accuracy is claimed here.
"""
from pathlib import Path
import argparse,json
from import_published import load
from projective_polynomial_inverse import bind
import dqi_budget_probe as probe
import prefix_high_copy_jacobi as prefix
import counter_quotient_jacobi as quotient
from cube_arithmetic import CubeArithmetic

HERE=Path(__file__).resolve().parent

def run(n=70,pairs=6,through=18,suffix=2,inverse_limit=1004):
 m=load();m.inv_word=bind(m,limit=inverse_limit)
 Base=m.Em
 class Counter(Base):
  def emit_mac(self,rt):self.b._cur_stats['toffoli']+=sum(t[0]=='P' for t in rt)
  def emit_free(self,rt):pass
 m.Em=Counter;b=m.Builder();b._items=probe.count_module().Builder()._items
 W=[b.qubits(12) for _ in range(n)]
 E=m.Em(b,*b.qubits(3),m.Alloc(b));math=CubeArithmetic(m,E)
 counts=prefix.schedule(n//2,pairs,through,suffix);stages=[]
 def stage(name,fn):
  before=b._cur_stats['toffoli'];fn()
  stages.append(dict(name=name,T=b._cur_stats['toffoli']-before,peak=b.peak_qubits))
 stage('extract',lambda:prefix.extract(m,E,math,W,counts=counts,through=through))
 with m.opened(E,'quotient_terminal',12) as c,m.opened(E,'quotient_metadata',(n//2).bit_length()+1) as metadata:
  stage('synthesize',lambda:quotient.synthesize(m,E,W,c,metadata))
  stage('synthesize_inverse',lambda:quotient.synthesize(m,E,W,c,metadata,True))
 stage('extract_inverse',lambda:prefix.extract(m,E,math,W,True,counts,through))
 core=b._cur_stats['toffoli'];scan=80*4095
 return dict(N=n,prefix=pairs,through=through,suffix=suffix,inverse_limit=inverse_limit,counts=counts,Q=b.peak_qubits,
             core_T=core,scan_T=scan,projected_total_T=core+scan,stages=stages,final_live=b.active_qubits,
             scope='exact emitted arithmetic counts; 80T scan projected, full native circuit unbuilt')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--pairs',type=int,nargs='+',default=[4,5,6,8]);p.add_argument('--through',type=int,nargs='+',default=[18]);p.add_argument('--suffix',type=int,nargs='+',default=[2]);p.add_argument('--inverse-limit',type=int,default=1004);p.add_argument('--out',default='prefix_counter_profiles');a=p.parse_args();rows=[]
 for b in a.pairs:
  for through in a.through:
   for suffix in a.suffix:
    r=run(pairs=b,through=through,suffix=suffix,inverse_limit=a.inverse_limit);rows.append(r);print(json.dumps(r),flush=True)
    (HERE/'results'/f'{a.out}.json').write_text(json.dumps(rows,indent=2)+'\n')
