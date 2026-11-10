"""Independent exact integer phase model; never imports a circuit builder."""
import hashlib, json, math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
MAX_X=3133452
TABLE=np.load(HERE/'reference/table.npy')

def rows(a,z,L,e):
    bit=lambda j: 0 if j<0 else (z>>min(j,L-1))&1
    for g in range((L+1)//2):
        digit=bit(2*g-1)+bit(2*g)-2*bit(2*g+1)
        value=np.zeros_like(z)
        for mag in (1,2):
            ee=e+2*g+(mag==2)
            rounded=a<<ee if ee>=0 else (a>>-ee)+((a>>(-ee-1))&1)
            value+=np.where(abs(digit)==mag,np.sign(digit)*rounded,0)
        yield value

def product(a,z,L,e):
    return sum(rows(a,z,L,e),np.zeros_like(z))

def evaluate(xs,detail=False):
    xs=np.asarray(xs,dtype=np.int64)
    # Integer bit-length, avoiding float log2 at exact powers of two.
    k=np.zeros_like(xs)
    for i in range(21):k+=xs<(1<<(21-i))
    m=xs<<k;s=k//2;addr=((m>>13)&255)+((k&1)<<8)
    c0,c1,c2=TABLE[addr].T
    z=(m&16383)^8192
    b=c1-product(c2,z,14,-11)
    rr=list(rows(b,z,14,-12));j=sum(rr,np.zeros_like(xs))
    p=c0-j;u=(1<<18)-j;a=c0-(1<<18)
    assert np.array_equal(p,a+u)
    assert np.all((u>=0)&(u<(1<<19)))
    assert np.array_equal(p,c0-sum(rr,np.zeros_like(xs)))
    result=(p<<s).astype(np.uint64)<<np.uint64(27)
    if detail:return result,dict(k=k,m=m,s=s,addr=addr,c0=c0,c1=c1,c2=c2,z=z,b=b,j=j,u=u,p=p)
    return result

def reference():
    path=HERE/'reference/exact_target.npz'
    if path.exists():
        with np.load(path) as d:return d['target'],d['tol']
    r=[math.isqrt((1<<130)//x) for x in range(1,MAX_X+1)]
    target=np.array([((v+1)>>1)&((1<<64)-1) for v in r],dtype=np.uint64)
    tol=np.array([v>>25 for v in r],dtype=np.uint64)
    np.savez(path,target=target,tol=tol)
    return target,tol

def synthesis(nrot,margin,tbase=0):
    eps=math.pi*margin/2**64/nrot if nrot else None
    return dict(rotation_count=nrot,rotation_eps=eps,
        synthesized_T=nrot*(.56*math.log2(1/eps)+5.3) if nrot else 0,
        synthesis_work=(tbase+(nrot*(.56*math.log2(1/eps)+5.3) if nrot else 0))/4)

def certify():
    xs=np.arange(1,MAX_X+1,dtype=np.int64);phase,d=evaluate(xs,True)
    target,tol=reference();used=np.abs((phase-target).view(np.int64)).astype(np.uint64)
    margin=int((tol.astype(np.int64)-used.astype(np.int64)).min());ratio=used/tol
    report=dict(status='PASS' if margin>0 else 'FAIL',scope='independent exact model; not native verification',
        inputs=MAX_X,worst_error_over_tolerance=float(ratio.max()),worst_x=int(ratio.argmax()+1),
        min_margin_units=margin,min_margin_turns=margin/2**64,
        ranges={name:[int(d[name].min()),int(d[name].max())] for name in ('b','j','u','p')},
        split_identity='C0-J=(C0-2^18)+(2^18-J)',
        table_sha256=hashlib.sha256((HERE/'reference/table.npy').read_bytes()).hexdigest())
    (HERE/'results/model.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    return report

if __name__=='__main__':certify()
