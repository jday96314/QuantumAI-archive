"""Exhaustive precision screen for the conditional fused-square build plan.

Reuse the accepted coefficient table, but evaluate linear term exactly and
round the coordinate only inside the quadratic correction. Omit/round square
columns before multiplying curvature. No circuit-resource claim is made.
"""
import json,math,time
from pathlib import Path
import numpy as np
from prototype import evaluate,TABLE,N
H=Path(__file__).resolve().parent

if __name__=='__main__':
    with np.load(H.parent/'inv-phase/reference/exact_target.npz') as d:target=d['target'];tol=d['tol']
    cases=[(q,b) for q in (0,2,3,4,5,6) for b in (0,10,12,14,16)]
    reports=[]
    for q,b in cases:
        minmargin=1<<64;maxratio=0;worst=0;maxabs=0
        for first in range(1,N+1,131072):
            xs=np.arange(first,min(N+1,first+131072),dtype=np.int64);k,m,s,p=evaluate(xs)
            c0,c1,c2=TABLE[((m>>13)&255)+((k&1)<<8)].T;z=(m&16383)-8192
            zz=((z+(1<<(q-1)))>>q)<<q if q else z
            maxabs=max(maxabs,int(np.abs(zz>>q).max()))
            square=zz*zz
            if b:square=((square+(1<<(b-1)))>>b)<<b
            num=(c0<<23)-(c1*z<<11)+c2*square
            phase=(num.astype(np.uint64)<<(s+4).astype(np.uint64))
            used=np.abs((phase-target[first-1:first-1+len(xs)]).view(np.int64)).astype(np.uint64)
            to=tol[first-1:first-1+len(xs)];ratio=used/to;j=int(ratio.argmax())
            if float(ratio[j])>maxratio:maxratio=float(ratio[j]);worst=int(xs[j])
            minmargin=min(minmargin,int((to.astype(np.int64)-used.astype(np.int64)).min()))
        report=dict(coordinate_low_bits_rounded=q,square_low_bits_rounded=b,
            signed_coordinate_bits=maxabs.bit_length()+1,unsigned_absolute_coordinate_bits=maxabs.bit_length(),
            maximum_absolute_coordinate=maxabs,square_significant_bits=27-max(2*q,b),
            inputs=N,status='MODEL_PASS' if minmargin>0 else 'MODEL_FAIL',
            worst_ratio=maxratio,worst_x=worst,min_margin_units=minmargin,min_margin_turns=minmargin/2**64)
        reports.append(report);print(json.dumps(report),flush=True)
    (H/'results/curvature_frontier.json').write_text(json.dumps(reports,indent=2)+'\n')
