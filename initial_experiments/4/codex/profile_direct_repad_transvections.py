"""Resource-only transvection profile for N70 squarefree cofactor padding."""
import json,time
from pathlib import Path
from import_published import load
import deferred_euclid_padding as padding
from controlled_linear import apply as linear


def main():
 rows=[]
 for fast in [False,True]:
  m=load();b=m.Builder();words=[b.qubits(12) for _ in range(70)];flags=b.qubits(35);E=m.Em(b,*b.qubits(3),m.Alloc(b));old=padding.emit_ops
  if fast:padding.emit_ops=lambda bb,ops,flag,inverse=False:linear(bb,ops,flag,inverse)
  start=time.time()
  try:padding.direct_repad(m,E,words,flags)
  finally:padding.emit_ops=old
  rows.append(dict(controlled_transvections=fast,Q=b.next_qubit,T=b._stream_totals()[0],ops=len(b._items),seconds=time.time()-start))
 report=dict(status='PROFILED',scope='resource-only common monomial to squarefree cofactor padding at N70; affine transvection helper independently tested, this full adapter is not native-verified',rows=rows)
 Path('results/direct_repad_transvection_profile.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
