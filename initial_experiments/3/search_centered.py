"""Search centered interpolation precision and fit each intercept interval."""
import json,time
from pathlib import Path
import numpy as np
import centered_invsqrt as p
HERE=Path(__file__).resolve().parent
original=p.table_for
rng=np.random.default_rng(101)
xs=np.unique(np.r_[np.arange(1<<20,1<<22,61),rng.integers(1<<20,1<<22,35000)])
records=[]
for f in [5,6,7,8]:
    addr=xs>>(20-f);starts=np.r_[0,np.flatnonzero(np.diff(addr))+1]
    for wc in [26,27,28,29]:
        ideal=2.**(wc+10)/np.sqrt(xs);tol=2.**-24
        for w1 in range(max(2,wc-f-5),wc-f+1):
            for w2 in range(max(2,wc-2*f-6),max(3,wc-2*f+1)):
                widths=(wc,w1,w2);v=p.model(xs,f,widths,0)
                lo=np.ceil(np.maximum.reduceat(ideal*(1-tol)-v,starts)).astype(int)
                hi=np.floor(np.minimum.reduceat(ideal*(1+tol)-v,starts)).astype(int)
                # Keep the largest normalized value below 1 to avoid a
                # carry beyond the intercept register's significant bits.
                hi[0]=min(hi[0],(1<<wc)-1-v[0])
                if (lo>hi).any():continue
                tab=original(f,widths,0).copy();tab[1<<f:,0]+=np.rint((lo+hi)/2).astype(int)
                def fit(ff,ww,bias=0):return tab if ff==f and tuple(ww)==widths else original(ff,ww,bias)
                p.table_for=fit
                try:c=p.build(f,widths,0)
                finally:p.table_for=original
                q=c.num_qubits;t=c.operation_statistics().num_toffoli_like
                rec=dict(frac=f,widths=widths,qubits=q,toffolis=t,estimated_score=q*(t+260),table=tab)
                records.append(rec)
                print('SAMPLE',json.dumps({k:v for k,v in rec.items() if k!='table'}),flush=True)
records.sort(key=lambda x:x['estimated_score'])
full=np.arange(1,3133453)
for idx,rec in enumerate(records[:20]):
    f=rec['frac'];widths=tuple(rec['widths']);tab=rec['table']
    def fit(ff,ww,bias=0):return tab if ff==f and tuple(ww)==widths else original(ff,ww,bias)
    p.table_for=fit
    try:
        vals=p.model(full,f,widths,0)
        err=vals/2.**(widths[0]+10)*np.sqrt(full)-1
        rec['max_relative_error']=float(abs(err).max())
        rec['full_pass']=rec['max_relative_error']<2.**-24
        if rec['full_pass']:
            path=HERE/f'candidates/centered_{idx}'
            c=p.build(f,widths,0);c.to_file(str(path.with_suffix('.kmb')),format='kmb');np.save(path.with_suffix('.npy'),tab)
            rec['circuit']=str(path.with_suffix('.kmb'));rec['table_path']=str(path.with_suffix('.npy'))
        print('FULL',json.dumps({k:v for k,v in rec.items() if k!='table'}),flush=True)
    finally:p.table_for=original
(HERE/'results/centered_search.json').write_text(json.dumps([{k:v for k,v in r.items() if k!='table'} for r in records],indent=2))
