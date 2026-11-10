"""All-domain phase quantization probes; no coefficient or segmentation search."""
import json
import numpy as np
from model import HERE,MAX_X,evaluate,rows,reference

def main():
    target,tol=reference()
    names=['exact_rows','exact_rational_quadratic']+[f'row_nearest_{j}' for j in range(1,6)]
    records={name:dict(model=name,max_error_ratio=0.,worst_x=0,min_margin_units=2**64) for name in names}
    for first in range(1,MAX_X+1,131072):
        xs=np.arange(first,min(first+131072,MAX_X+1),dtype=np.int64)
        exact,d=evaluate(xs,True);s=d['s'];z=np.where(d['z']>=8192,d['z']-16384,d['z'])
        # Replace BOTH rounded products by exact rational products; it is a
        # different approximation and is not assumed to preserve acceptance.
        numerator=(d['c0']<<23)-(d['c1']*z<<11)+d['c2']*z*z
        phases={'exact_rows':exact,'exact_rational_quadratic':numerator.astype(np.uint64)<<(s+4).astype(np.uint64)}
        rr=list(rows(d['b'],d['z'],14,-12))
        for bits in range(1,6):
            unit=1<<bits
            rounded=[np.sign(v)*((np.abs(v)+(unit//2))//unit)*unit for v in rr]
            p=d['c0']-sum(rounded,np.zeros_like(xs))
            phases[f'row_nearest_{bits}']=(p<<s).astype(np.uint64)<<np.uint64(27)
        for name,phase in phases.items():
            used=np.abs((phase-target[first-1:first-1+len(xs)]).view(np.int64)).astype(np.uint64)
            allowed=tol[first-1:first-1+len(xs)];ratio=used/allowed;j=int(ratio.argmax());r=records[name]
            if ratio[j]>r['max_error_ratio']:r.update(max_error_ratio=float(ratio[j]),worst_x=int(xs[j]))
            r['min_margin_units']=min(r['min_margin_units'],int((allowed.astype(np.int64)-used.astype(np.int64)).min()))
    out=[]
    for name in names:
        r=records[name];r.update(status='MODEL_PASS' if r['min_margin_units']>0 else 'MODEL_FAIL',inputs=MAX_X)
        out.append(r);print(json.dumps(r))
    (HERE/'results/error_frontier.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
