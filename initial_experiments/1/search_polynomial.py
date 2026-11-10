"""Search interpolation width and precision, then exhaustively check finalists."""
import itertools,json,time
from pathlib import Path
import numpy as np
from polynomial_invsqrt import model,build
rng=np.random.default_rng(20261023)
xs=np.unique(np.r_[np.arange(1,4097),rng.integers(1,3133453,80000),np.arange(1048576,3133453,64)])
s=(22-(np.floor(np.log2(xs)).astype(np.int64)+1))//2
root=np.sqrt(xs)/2.**(10-s)
target=2.**-24
records=[]
for order,frac in [(2,5),(2,6),(2,7),(3,3),(3,4),(3,5),(4,2),(4,3)]:
 for wc0 in range(27,31):
  for extra in [-1,0,1,2]:
   widths=tuple([wc0]+[max(2,wc0-j*frac+extra) for j in range(1,order+1)])
   try: vals=model(xs,frac,widths,0)>>s
   except (ValueError,AssertionError):continue
   exact=2.**wc0/root
   lo=int(np.ceil(np.max(exact*(1-target)-vals)))
   hi=int(np.floor(np.min(exact*(1+target)-vals)))
   if lo>hi:continue
   bias=round((lo+hi)/2)
   c=build(frac,widths,bias)
   toff=c.operation_statistics().num_toffoli_like
   est=c.num_qubits*(toff+260)
   rec=dict(frac=frac,widths=widths,bias=bias,qubits=c.num_qubits,toffolis=toff,estimated_score=est)
   records.append(rec);print(json.dumps(rec),flush=True)
records.sort(key=lambda d:d['estimated_score'])
Path('initial_experiments/results/polynomial-search.json').write_text(json.dumps(records,indent=2))
fullxs=np.arange(1,3133453)
for i,r in enumerate(records[:12]):
 vals=model(fullxs,r['frac'],r['widths'],r['bias'])
 err=vals/2.**(r['widths'][0]+10)*np.sqrt(fullxs)-1
 r['full_max_relative_error']=float(abs(err).max())
 r['full_pass']=r['full_max_relative_error']<target
 print('FULL',json.dumps(r),flush=True)
 if r['full_pass']:
  path=Path(f'initial_experiments/candidates/poly-search-{i}.kmb')
  build(r['frac'],r['widths'],r['bias']).to_file(path,format='kmb');r['circuit']=str(path)
Path('initial_experiments/results/polynomial-search.json').write_text(json.dumps(records,indent=2))
