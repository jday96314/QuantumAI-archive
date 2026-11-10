"""Fit one constant correction per interpolation interval, exhaustively."""
import json
from pathlib import Path
import numpy as np
import polynomial_invsqrt as p
original=p.table_for
xs=np.arange(1<<20,1<<22,dtype=np.int64)
results=[]
for order,frac in [(2,5),(2,6),(2,7),(3,3),(3,4),(3,5),(4,2),(4,3)]:
 for wc in range(26,30):
  for extra in [-2,-1,0,1]:
   widths=tuple([wc]+[max(2,wc-j*frac+extra) for j in range(1,order+1)])
   try:vals=p.model(xs,frac,widths,0)
   except (ValueError,AssertionError):continue
   ideal=2.**(wc+10)/np.sqrt(xs)
   step=1<<(20-frac);target=2.**-24
   lo=np.ceil((ideal*(1-target)-vals).reshape(-1,step).max(axis=1)).astype(int)
   hi=np.floor((ideal*(1+target)-vals).reshape(-1,step).min(axis=1)).astype(int)
   if (lo>hi).any():continue
   biases=np.rint((lo+hi)/2).astype(int)
   tab=original(frac,widths,0).copy();tab[1<<frac:,0]+=biases
   if (tab<0).any():continue
   def corrected(f,w,b=0):
    if f==frac and tuple(w)==widths:return tab
    return original(f,w,b)
   p.table_for=corrected
   try:
    c=p.build(frac,widths,0);toff=c.operation_statistics().num_toffoli_like
    name=f'poly-interval-f{frac}-w'+ '-'.join(map(str,widths))
    record=dict(frac=frac,widths=widths,qubits=c.num_qubits,toffolis=toff,estimated_score=c.num_qubits*(toff+250),table=f'initial_experiments/candidates/{name}.npy',circuit=f'initial_experiments/candidates/{name}.kmb')
    np.save(record['table'],tab);c.to_file(record['circuit'],format='kmb')
    results.append(record);print(json.dumps(record),flush=True)
   finally:p.table_for=original
results.sort(key=lambda r:r['estimated_score'])
Path('initial_experiments/results/interval-bias-search.json').write_text(json.dumps(results,indent=2))
