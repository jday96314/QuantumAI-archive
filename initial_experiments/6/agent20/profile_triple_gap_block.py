"""Actual single-block emitter counts; no full decoder claim."""
import json
from pathlib import Path
from import_published import load
from projective_polynomial_inverse import bind
import dqi_budget_probe as probe
import conditional_ccz_builder
from compact_blocked_cfraction import Arithmetic
import triple_gap_block as triple
HERE=Path(__file__).resolve().parent
rows=[]
for limit in (935,1004):
 m=load();conditional_ccz_builder.install(m);m.inv_word=bind(m,limit=limit);Base=m.Em
 class Counter(Base):
  def emit_mac(self,rt):self.b._cur_stats['toffoli']+=sum(t[0]=='P' for t in rt)
  def emit_free(self,rt):pass
 m.Em=Counter;b=m.Builder();b._items=probe.count_module().Builder()._items;W=[b.qubits(12) for _ in range(70)];E=m.Em(b,*b.qubits(3),m.Alloc(b));triple.block(m,E,Arithmetic(m,E),W)
 r=dict(N=70,inverse_limit=limit,Q=b.peak_qubits,T=b._cur_stats['toffoli'],scope='single shared triple-gap block, no complete extraction/synthesis/scan');rows.append(r);print(json.dumps(r),flush=True);(HERE/'results/triple_gap_block_profiles.json').write_text(json.dumps(rows,indent=2)+'\n')
