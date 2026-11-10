"""Exact emitted block profiles; no complete decoder or native acceptance."""
import json,time
from pathlib import Path
import dqi_budget_probe as probe
import general_jacobi as general
import continuant_head_frame as compact
from import_published import load
from projective_polynomial_inverse import bind
from cube_arithmetic import CubeArithmetic
HERE=Path(__file__).resolve().parent


def run(N,b,kind):
 m=load();m.inv_word=bind(m)
 Base=m.Em
 class Counter(Base):
  def emit_mac(self,rt):self.b._cur_stats['toffoli']+=sum(t[0]=='P' for t in rt)
  def emit_free(self,rt):pass
 m.Em=Counter;bb=m.Builder();bb._items=probe.count_module().Builder()._items;W=[bb.qubits(12) for _ in range(N)]
 E=m.Em(bb,*bb.qubits(3),m.Alloc(bb));math=CubeArithmetic(m,E)
 old=general.parameters
 if kind=='compact':general.parameters=compact.parameters
 before=bb._cur_stats['toffoli'];general.block(m,E,math,W,b);T=bb._cur_stats['toffoli']-before;general.parameters=old
 return dict(N=N,pairs=b,kind=kind,Q=bb.peak_qubits,T=T,scope='emitted single block, forward only, no tape retained')
if __name__=='__main__':
 rows=[]
 for b in (2,3,4):
  for kind in ('prior','compact'):
   r=run(70,b,kind);rows.append(r);print(json.dumps(r),flush=True);(HERE/'results/continuant_frame_direct_block_profiles.json').write_text(json.dumps(rows,indent=2)+'\n')
