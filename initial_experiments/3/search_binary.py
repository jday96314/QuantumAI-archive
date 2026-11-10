"""Fit binary-normalized interpolation tables, then check every legal input."""
import json
from pathlib import Path
import numpy as np
import binary_invsqrt as p
HERE=Path(__file__).resolve().parent
original=p.table_for;rng=np.random.default_rng(103)
# Both normalization parity classes, sampled in normalized mantissa order.
norm=np.unique(np.r_[np.arange(1<<21,1<<22,53),rng.integers(1<<21,1<<22,25000)])
xs=np.r_[norm,norm//2]
records=[]
for f in [5,6,7,8]:
    k=21-np.floor(np.log2(xs)).astype(np.int64)
    addr=(((xs<<k)>>(21-f))&((1<<f)-1))+((k&1)<<f)
    order=np.argsort(addr);xx=xs[order];addr=addr[order]
    starts=np.r_[0,np.flatnonzero(np.diff(addr))+1]
    for wc in [26,27,28,29]:
        ideal=2.**(wc+10)/np.sqrt(xx);ideal/=2.**(k[order]//2)
        for w1 in range(max(2,wc-f-5),wc-f+2):
            for w2 in range(max(2,wc-2*f-6),max(3,wc-2*f+3)):
                widths=(wc,w1,w2)
                vv=p.model(xx,f,widths,0)//(1<<(k[order]//2))
                tol=2.**-24
                lo=np.ceil(np.maximum.reduceat(ideal*(1-tol)-vv,starts)).astype(int)
                hi=np.floor(np.minimum.reduceat(ideal*(1+tol)-vv,starts)).astype(int)
                if (lo>hi).any():continue
                tab=original(f,widths,0).copy();tab[:,0]+=np.rint((lo+hi)/2).astype(int)
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
        rec['max_relative_error']=float(abs(err).max());rec['full_pass']=rec['max_relative_error']<2.**-24
        if rec['full_pass']:
            path=HERE/f'candidates/binary_{idx}'
            c=p.build(f,widths,0);c.to_file(str(path.with_suffix('.kmb')),format='kmb');np.save(path.with_suffix('.npy'),tab)
            rec.update(circuit=str(path.with_suffix('.kmb')),table_path=str(path.with_suffix('.npy')))
        print('FULL',json.dumps({k:v for k,v in rec.items() if k!='table'}),flush=True)
    finally:p.table_for=original
(HERE/'results/binary_search.json').write_text(json.dumps([{k:v for k,v in r.items() if k!='table'} for r in records],indent=2))
