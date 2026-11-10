"""Exact current polynomial reconstruction resources in two live layouts."""
from pathlib import Path
import json
from import_published import load
from projective_polynomial_inverse import bind
import generalized_gap_quotient as general
import dqi_budget_probe as probe
HERE=Path(__file__).resolve().parent
rows=[]
for layout,n,raw,metadata,positions in [('full_triples',35,0,7,[i for i in range(34) if i%3!=2]),('six_head',6,72,0,[0,2,4])]:
 for limit in (935,1004,1024):
  m=load();m.inv_word=bind(m,limit=limit);Base=m.Em
  class Counter(Base):
   def emit_mac(self,rt):self.b._cur_stats['toffoli']+=sum(t[0]=='P' for t in rt)
   def emit_free(self,rt):pass
  m.Em=Counter;b=m.Builder();b._items=probe.count_module().Builder()._items;W=[b.qubits(12) for _ in range(2*n)]
  if layout=='six_head':b.qubits(840-24*n)
  if raw:b.qubits(raw)
  c=b.qubits(12);bits=b.qubits(len(positions));lookup=dict(zip(positions,bits));gaps=[lookup.get(i) for i in range(n)]
  if metadata:b.qubits(metadata)
  E=m.Em(b,*b.qubits(3),m.Alloc(b));F,G=general.base.layout(W,c);general.synthesize(m,E,F,G,gaps)
  r=dict(layout=layout,n=n,inverse_limit=limit,Q=b.peak_qubits,T=b._cur_stats['toffoli'],gap_flags=len(bits),metadata_bits=metadata,scope='polynomial reconstruction only with explicit live data/metadata/prefix-copy footprint; no Jacobi conversion, normalization or full decoder');rows.append(r);print(json.dumps(r),flush=True);(HERE/'results/general_gap_quotient3_profiles.json').write_text(json.dumps(rows,indent=2)+'\n')
