"""Complete-domain exact fixed-point phase certificate, independent of kickmix.

Integer Horner follows each individual signed Booth row, including its rounding.
Reference phase/tolerance use the host's isqrt formula, rather than float error.
"""
import argparse, collections, json, math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
MAX_X=3133452

def product(a,z,L,e):
    out=np.zeros_like(z)
    bit=lambda k:0 if k<0 else (z>>min(k,L-1))&1
    for g in range((L+1)//2):
        digit=bit(2*g-1)+bit(2*g)-2*bit(2*g+1)
        for mag in (1,2):
            ee=e+2*g+(mag==2)
            v=a<<ee if ee>=0 else (a>>-ee)+((a>>(-ee-1))&1)
            out+=np.where(abs(digit)==mag,np.sign(digit)*v,0)
    return out

def evaluate(model,xs,states=False):
    wc,w1,w2=model['widths'];xs=np.asarray(xs,dtype=np.int64)
    k=21-np.floor(np.log2(xs)).astype(np.int64);n=xs<<k
    y=np.zeros_like(xs);slope_range=[2**63,-1];y_range=[2**63,-1]
    max_products=[0,0];sign_violations=[0,0]
    depth_count=collections.Counter()
    for p,lo,d,f,degree,a,b,c,err in model['leaves']:
        mask=((k&1)==p)&(n>=2**21+lo)&(n<2**21+lo+2**(21-d))
        if not mask.any():continue
        z=(n[mask]&((1<<(21-f))-1))^(1<<(20-f))
        first=product(c,z,21-f,w1-w2-21) if degree==2 else np.zeros_like(z)
        slope=b-first
        second=product(slope,z,21-f,wc-w1-21);value=a-second
        if states:
            signed=np.where(z>=(1<<(20-f)),z-(1<<(21-f)),z)
            for j,pp in enumerate([second,first]):
                max_products[j]=max(max_products[j],int(np.abs(pp).max()))
                sign_violations[j]+=int(((pp*signed)<0).sum())
        slope_range=[min(slope_range[0],int(slope.min())),max(slope_range[1],int(slope.max()))]
        y_range=[min(y_range[0],int(value.min())),max(y_range[1],int(value.max()))]
        y[mask]=value;depth_count[d]+=int(mask.sum())
    result=y<<(k//2)
    if states:
        base=min(l[2] for l in model['leaves']);split=[wc-base-2,w1-base-2]
        return result,dict(slope_range=slope_range,mantissa_range=y_range,depth_inputs=dict(depth_count),
            max_product_magnitudes=max_products,product_sign_violations=sign_violations,
            signed_accumulator_splits=split,
            signed_accumulator_magnitudes_fit=all(v<(1<<s) for v,s in zip(max_products,split)))
    return result

def reference():
    path=HERE/'reference/target_tol.npz'
    if path.exists():
        with np.load(path) as d:return d['target'],d['tol']
    rs=[math.isqrt((1<<130)//x) for x in range(1,MAX_X+1)]
    target=np.array([((r+1)>>1)&((1<<64)-1) for r in rs],dtype=np.uint64)
    tol=np.array([r>>25 for r in rs],dtype=np.uint64)
    np.savez(path,target=target,tol=tol)
    return target,tol

def certificate(model):
    xs=np.arange(1,MAX_X+1,dtype=np.int64);wc=model['widths'][0]
    y,states=evaluate(model,xs,True)
    k=21-np.floor(np.log2(xs)).astype(np.int64);n=xs<<k
    base=min(l[2] for l in model['leaves'])
    residual=(y>>(k//2))+(((n>>(21-base))&((1<<base)-1))<<(wc-base-2))-(11<<(wc-4))-((k&1)<<(wc-2))
    states['affine_residual_range']=[int(residual.min()),int(residual.max())]
    states['affine_base']=base
    phase=y.astype(np.uint64)<<np.uint64(64-wc-10)
    target,tol=reference();error=(phase-target).view(np.int64)
    used=np.abs(error).astype(np.uint64);valid=used<tol
    ratio=used/tol;idx=int(ratio.argmax())
    margin=int((tol.astype(np.int64)-used.astype(np.int64)).min())
    rotations=wc+7;eps=math.pi*margin/2**64/rotations
    synth=(1+rotations*(.56*math.log2(1/eps)+5.3))/4 if eps>0 else None
    return dict(status='PASS' if valid.all() else 'FAIL',scope='exact integer model, not native verification',
                num_inputs=MAX_X,worst_x=idx+1,worst_error_over_tolerance=float(ratio[idx]),
                min_margin_units=margin,min_margin_turns=margin/2**64,rotation_count=rotations,
                rotation_eps=eps,synthesis_toffoli_equivalent=synth,
                leaf_count=len(model['leaves']),depth_histogram=dict(collections.Counter(l[2] for l in model['leaves'])),
                degree_histogram=dict(collections.Counter(l[4] for l in model['leaves'])),**states)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('models',nargs='+',type=Path);a=p.parse_args()
    for path in a.models:
        report=certificate(json.loads(path.read_text()));report['model']=str(path)
        (HERE/'results'/f'{path.stem}_certificate.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report),flush=True)
