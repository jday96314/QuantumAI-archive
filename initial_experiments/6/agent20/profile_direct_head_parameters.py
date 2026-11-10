"""Bounded cache/width tradeoff for compact head affine preparation."""
import json
from pathlib import Path
from import_published import load
from projective_polynomial_inverse import bind
from cube_arithmetic import CubeArithmetic
from field_dag import TreeEvaluator,Evaluator
import continuant_head_frame as compact
import dqi_budget_probe as probe
HERE=Path(__file__).resolve().parent

def run(b,capacity):
 m=load();m.inv_word=bind(m);Base=m.Em
 class Counter(Base):
  def emit_mac(self,rt):self.b._cur_stats['toffoli']+=sum(t[0]=='P' for t in rt)
  def emit_free(self,rt):pass
 m.Em=Counter;bb=m.Builder();bb._items=probe.count_module().Builder()._items;W=[bb.qubits(12) for _ in range(70)];E=m.Em(bb,*bb.qubits(3),m.Alloc(bb));math=CubeArithmetic(m,E)
 class Cached(Evaluator):
  def __init__(self,*args,**kwargs):super().__init__(*args,capacity=capacity,**kwargs)
 with compact.parameters(m,E,math,W[:2*b],evaluator=TreeEvaluator if capacity is None else Cached):pass
 return dict(pairs=b,capacity=capacity,Q=bb.peak_qubits,T=bb._cur_stats['toffoli'],scope='full parameter preparation/erasure/frame at840 input qubits; count emitter')
if __name__=='__main__':
 rows=[]
 for b in (2,3,4):
  for cap in (None,4,6,8,10,16):
   r=run(b,cap);rows.append(r);print(json.dumps(r),flush=True);(HERE/'results/direct_head_parameter_profiles.json').write_text(json.dumps(rows,indent=2)+'\n')
