"""Search independent slope/curvature widths with per-interval integer correction."""
import json
from pathlib import Path
import numpy as np
import polynomial_invsqrt as p
original=p.table_for
rng=np.random.default_rng(20261025)
xs=np.unique(np.r_[np.arange(1<<20,1<<22,37),rng.integers(1<<20,1<<22,70000)])
results=[]
for frac in [6,7,8]:
    addresses=xs>>(20-frac);starts=np.r_[0,np.flatnonzero(np.diff(addresses))+1]
    idealbase=1024/np.sqrt(xs)
    for wc in range(26,30):
        ideal=2.**wc*idealbase
        for w1 in range(wc-frac-3,wc-frac+2):
            for w2 in range(max(2,wc-2*frac-2),wc-2*frac+3):
                widths=(wc,w1,w2)
                vals=p.model(xs,frac,widths,0)
                lo=np.ceil(np.maximum.reduceat(ideal*(1-2.**-24)-vals,starts)).astype(int)
                hi=np.floor(np.minimum.reduceat(ideal*(1+2.**-24)-vals,starts)).astype(int)
                if (lo>hi).any():continue
                tab=original(frac,widths,0).copy();tab[1<<frac:,0]+=np.rint((lo+hi)/2).astype(int)
                def corrected(f,w,b=0):return tab if f==frac and tuple(w)==widths else original(f,w,b)
                p.table_for=corrected
                try:c=p.build(frac,widths,0)
                finally:p.table_for=original
                t=c.operation_statistics().num_toffoli_like
                rec=dict(frac=frac,widths=widths,qubits=c.num_qubits,toffolis=t,estimated_score=c.num_qubits*(t+250),table_data=tab.tolist())
                results.append(rec);print('SAMPLE',json.dumps({k:v for k,v in rec.items() if k!='table_data'}),flush=True)
results.sort(key=lambda r:r['estimated_score'])
fullxs=np.arange(1<<20,1<<22)
for i,rec in enumerate(results[:25]):
    frac=rec['frac'];widths=tuple(rec['widths']);tab=np.array(rec['table_data'],dtype=np.int64)
    def corrected(f,w,b=0):return tab if f==frac and tuple(w)==widths else original(f,w,b)
    p.table_for=corrected
    try:
        vals=p.model(fullxs,frac,widths,0)
        err=vals/2.**(widths[0]+10)*np.sqrt(fullxs)-1
        rec['full_max_relative_error']=float(abs(err).max())
        rec['full_pass']=rec['full_max_relative_error']<2.**-24
        if rec['full_pass']:
            path=Path(f'initial_experiments/candidates/poly-mixed-{i}')
            c=p.build(frac,widths,0);c.to_file(path.with_suffix('.kmb'),format='kmb');np.save(path.with_suffix('.npy'),tab)
            rec['circuit']=str(path.with_suffix('.kmb'));rec['table']=str(path.with_suffix('.npy'))
        print('FULL',json.dumps({k:v for k,v in rec.items() if k!='table_data'}),flush=True)
    finally:p.table_for=original
Path('initial_experiments/results/mixed-precision-search.json').write_text(json.dumps([{k:v for k,v in r.items() if k!='table_data'} for r in results],indent=2))
