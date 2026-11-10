"""Reproduce compiled padded component widths, without claiming a full decoder."""
import json
from pathlib import Path
from local_schur import emit as head
from routed_cofactor import emit as cofactor
from shared_slot import emit as slot
from moment_probe import emit as moment
rows=[]
for name,make in [('routed_schur',lambda:head('forward',True)),('shared_cofactor',lambda:cofactor(19,'forward',3,True)),('carry_finish_slot',lambda:slot(18,'forward',True)),('conditional_moment9',lambda:moment(19,9,'forward',True))]:
    c=make();rows.append(dict(name=name,Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like))
report=dict(scope='Compiled native arithmetic embeddings with live operand padding; not a full decoder; pending-prefix placement and controller construction excluded',rows=rows,max_Q=max(r['Q'] for r in rows))
Path('results/embedded_components.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
