"""Emitted resource frontier for the quotient-frame prefix blocks."""
from pathlib import Path
import json
from import_published import load
from projective_polynomial_inverse import bind
import dqi_budget_probe as probe
from compact_blocked_cfraction import Arithmetic
import prefix_lookahead_quotient_block as block
HERE=Path(__file__).resolve().parent
rows=[]
for pairs in (4,6,8,10):
 for limit in (1004,1024):
  m=load();m.inv_word=bind(m,limit=limit);Base=m.Em
  class Counter(Base):
   def emit_mac(self,rt):self.b._cur_stats['toffoli']+=sum(t[0]=='P' for t in rt)
   def emit_free(self,rt):pass
  m.Em=Counter;b=m.Builder();b._items=probe.count_module().Builder()._items;W=[b.qubits(12) for _ in range(70)];E=m.Em(b,*b.qubits(3),m.Alloc(b));block.block(m,E,Arithmetic(m,E),W,pairs)
  r=dict(N=70,pairs=pairs,inverse_limit=limit,Q=b.peak_qubits,T=b._cur_stats['toffoli'],scope='one prefix block with cofactor gauge retained and raw-high boundary reused; no full decoder');rows.append(r);print(json.dumps(r),flush=True);(HERE/'results/prefix_lookahead_quotient_profiles.json').write_text(json.dumps(rows,indent=2)+'\n')
