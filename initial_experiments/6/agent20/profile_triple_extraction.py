"""Exact whole extraction resources with three-slot gap blocks."""
from pathlib import Path
import json
from import_published import load
from projective_polynomial_inverse import bind
import dqi_budget_probe as probe
from compact_blocked_cfraction import Arithmetic
import triple_gap_block as triple
import lookahead_jacobi as look
import compact_lookahead_parameters as compact
HERE=Path(__file__).resolve().parent
m=load();m.inv_word=bind(m,limit=935);Base=m.Em
class Counter(Base):
 def emit_mac(self,rt):self.b._cur_stats['toffoli']+=sum(t[0]=='P' for t in rt)
 def emit_free(self,rt):pass
m.Em=Counter;b=m.Builder();b._items=probe.count_module().Builder()._items;W=[b.qubits(12) for _ in range(70)];E=m.Em(b,*b.qubits(3),m.Alloc(b));math=Arithmetic(m,E);look.parameters=compact.parameters;rows=[]
for k in range(0,35,3):
 words=W[2*k:];before=b._cur_stats['toffoli']
 if len(words)==2:look.jacobi.pair(m,E,math,words)
 elif len(words)==4:look.block(m,E,math,words)
 else:triple.block(m,E,math,words)
 rows.append(dict(start=k,T=b._cur_stats['toffoli']-before,peak=b.peak_qubits));print(json.dumps(rows[-1]),flush=True)
report=dict(N=70,Q=b.peak_qubits,extract_T=b._cur_stats['toffoli'],stages=rows,scope='complete triple extraction only; quotient conversion/reconstruction, inverse extraction and scan omitted')
(HERE/'results/triple_extraction_profile.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
