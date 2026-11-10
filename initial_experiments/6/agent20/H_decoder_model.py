"""Independent law-sampled H endpoint check using the compact reference frames.

The reference decoder runs only on syndromes. Constructed errors are used to
check locator/evaluator phases, not to construct the candidate transform.
Read-only serial carried-gap model stands in for equivalent hybrid arithmetic;
native complete integration remains a separate required check.
"""
import os,sys,json,time
from pathlib import Path
import numpy as np
from import_published import load
from sylvester_padding_screen import sylvester,prod
from dqi.verify_submission_dqi import optimal_dqi_weights
from published_dqi import isomorphism

HERE=Path(__file__).resolve().parent;GAP=HERE/'../../5/dqi-gap'
sys.path.insert(0,str(GAP.resolve()));oldcwd=Path.cwd();os.chdir(GAP)
try:import serial_reference as reference
finally:os.chdir(oldcwd)

def evaluate(m,p,x):
    value=0
    for coef in reversed(p):value=m.gmul(value,x)^coef
    return value

def wronskian(m,a,b,padding,repair=True):
    sq=lambda x:m.gpow(x,2048);n=len(a);t=n//2
    C=list(map(sq,a[::2]));D=list(map(sq,a[1::2]));E=list(map(sq,b[::2]));O=list(map(sq,b[1::2]))
    state=[False,padding//2] if repair else None
    h=sylvester(E,O,D,C[:-1],state)+[C[-1]]
    for i,x in enumerate(O[:-1]):h[t+i]^=m.gmul(C[-1],x)
    direct=[x^y for x,y in zip(prod(C,O,n),prod(D,E,n))]
    return h,direct,bool(state and state[0])

def run(shots=2000):
    m=load();rng=np.random.default_rng(2026110701);weights=optimal_dqi_weights(4095,35,12,2016)**2
    fw,bk,tr=isomorphism(m)
    assert all(tr(reference.mul(a,b),fw)==m.gmul(tr(a,fw),tr(b,fw)) for a,b in rng.integers(0,4096,(100,2)))
    started=time.time();counts=dict(base_wrong=0,degree_mismatch=0,high_A_nonzero=0,H_wrong=0,phase_wrong=0,repaired=0)
    failures=[];ds={}
    for shot in range(shots):
        d=int(rng.choice(36,p=weights));ds[d]=ds.get(d,0)+1
        points=rng.choice(np.arange(1,4096),d,replace=False).tolist();values=rng.integers(1,4096,d).tolist()
        synd=[0]*70
        for a,e in zip(points,values):
            power=a
            for i in range(70):synd[i]^=m.gmul(e,power);power=m.gmul(power,a)
        bp,b,p,meta=reference.decode([tr(x,bk) for x in synd])
        bp=[tr(x,fw) for x in bp];b=[tr(x,fw) for x in b];p=tr(p,fw)
        pad=0
        for k,flag in enumerate(meta['flags']):
            if flag:pad^=(35-k)^(34-k)
        A=[m.gmul(p,x) for x in bp]
        degree_ok=(pad==35-d and not any(b[d+1:]) and b[d]==1)
        high_ok=not any(A[d:]);counts['degree_mismatch']+=not degree_ok;counts['high_A_nonzero']+=not high_ok
        rawA=A[-pad:]+A[:-pad] if pad else A
        # Explicit leading coefficient rotates with the full B bank.
        rawB=b[-pad:]+b[:-pad] if pad else b
        h,direct,used=wronskian(m,rawA,rawB,pad);counts['repaired']+=used
        hh=h[pad:]+h[:pad] if pad else h
        hbad=h!=direct;counts['H_wrong']+=hbad;basebad=False;phasebad=False
        for a,e in zip(points,values):
            beta=m.ginv(a);der=0
            for i in range(1,len(b),2):der^=m.gmul(b[i],m.gpow(beta,i-1))
            denom=m.gmul(evaluate(m,A,beta),der)
            base=m.gmul(m.gpow(beta,70),m.ginv(denom))
            hv=evaluate(m,hh,beta);new=m.gmul(m.gpow(beta,70),m.ginv(m.gmul(hv,hv)))
            basebad|=base!=e or evaluate(m,b,beta)!=0
            phasebad|=new!=e or evaluate(m,b,beta)!=0
        counts['base_wrong']+=basebad;counts['phase_wrong']+=phasebad
        if hbad or basebad or phasebad or not degree_ok or not high_ok:
            failures.append(dict(shot=shot,degree=d,padding=pad,H_wrong=hbad,base_wrong=basebad,phase_wrong=phasebad,meta=meta))
        if shot%250==0:print(json.dumps(dict(shot=shot,counts=counts,seconds=time.time()-started)),flush=True)
    report=dict(status='PASS_MODEL' if not any(v for k,v in counts.items() if k!='repaired') else 'MODEL_RESIDUAL_FAILURES',
                scope='native-field independent law stream; serial carried-gap decoder reference, raw H prepare and unpad, all physical root magnitudes',
                shots=shots,seed=2026110701,counts=counts,degrees=ds,failures=failures,seconds=time.time()-started)
    (HERE/'results/H_decoder_model.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':run()
