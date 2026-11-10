"""Exhaustive interval fitting of nearby integer slope/curvature coefficients."""
import json
from pathlib import Path
import numpy as np
import polynomial_invsqrt as p
original=p.table_for
frac=7;step=1<<(20-frac)
delta=np.arange(step,dtype=np.int64)
bitrows=[(delta>>(20-m))&1 for m in range(frac+1,21)]
results=[]
for widths in [(27,19,11),(27,20,12),(27,20,13),(28,19,11),(28,19,12),(28,20,11),(28,20,12)]:
    wc,w1,w2=widths;tab=original(frac,widths,0).copy();worst=0.
    for address in range(1<<frac,1<<(frac+2)):
        ideal=2.**(wc+10)/np.sqrt((address<<(20-frac))+delta)
        best=np.inf;chosen=None
        for ds in range(-3,4):
            c1=int(tab[address,1])+ds
            for dc in range(-3,4):
                c2=int(tab[address,2])+dc
                if c1<0 or c2<0:continue
                slope=np.full(step,c1,dtype=np.int64)
                for m,bit in zip(range(frac+1,21),bitrows):
                    e=w1-w2-m
                    slope-=bit*(c2<<e if e>=0 else c2>>-e)
                val=np.full(step,int(tab[address,0]),dtype=np.int64)
                for m,bit in zip(range(frac+1,21),bitrows):
                    e=wc-w1-m
                    val-=bit*(slope<<e if e>=0 else slope>>-e)
                diff=ideal-val
                center=int(round((float(diff.max())+float(diff.min()))/2))
                for bias in [center-1,center,center+1]:
                    error=float(np.max(np.abs((val+bias)/ideal-1)))
                    if error<best:best=error;chosen=(int(tab[address,0])+bias,c1,c2)
        tab[address]=chosen;worst=max(worst,best)
    print('FIT',widths,'max_error',worst,flush=True)
    if worst>=2.**-24:continue
    def corrected(f,w,b=0):return tab if f==frac and tuple(w)==widths else original(f,w,b)
    p.table_for=corrected
    try:c=p.build(frac,widths,0)
    finally:p.table_for=original
    name='poly-integer-'+'-'.join(map(str,widths));path=Path('initial_experiments/candidates')/name
    np.save(path.with_suffix('.npy'),tab);c.to_file(path.with_suffix('.kmb'),format='kmb')
    rec=dict(frac=frac,widths=widths,qubits=c.num_qubits,toffolis=c.operation_statistics().num_toffoli_like,max_relative_error=worst,circuit=str(path.with_suffix('.kmb')),table=str(path.with_suffix('.npy')))
    print('PASS',json.dumps(rec),flush=True);results.append(rec)
Path('initial_experiments/results/integer-polynomial-search.json').write_text(json.dumps(results,indent=2))
