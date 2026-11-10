"""Native-field classical reference for the carried-gap compact representation."""
from classical_repair import mul,inv,EX,trim,divide,product,plus,monic
from check_gap_serial import reference


def decode(s):
    n=len(s);w=n//2;s=list(s);waits=[];flags=[];terminal=0
    wait=int(s[-1]==0 and s[-2]!=0);waits.append(wait)
    if wait:bp=[1];b=[0,1];c=s[:-2];p=1;stash=[s[-1],s[-2]];flags.append(0)
    else:
        flag=int(s[-1]==0);flags.append(flag);s[-1]^=flag
        p=inv(s[-1]);al=mul(s[-2],p);bp=[1];b=[al,1];c=[x^mul(al,y) for x,y in zip(s[:-2],s[1:-1])];stash=[0,0]
    for k in range(1,w):
        carry=waits[-1];shift=(waits[k-2] if k>=2 else 0)^terminal
        wait=int(not(c[-1] or not c[-2] or carry or shift))
        bp,b,c,p,stash,flag=reference(bp,b,c,p,stash,carry,wait,shift)
        terminal^=(waits[k-2] if k>=2 else 0)&flag
        waits.append(wait);flags.append(flag)
    shift=(waits[-2] if w>=2 else 0)^terminal
    count=0
    for k,f in enumerate(flags):
        if f:count^=(w-k)^(w-k-1)
    def rotate(v,a):
        a%=len(v);return v[a:]+v[:a]
    b=rotate(b[:-1],count)+[1];bp=rotate(bp,count)
    for k,f in enumerate(flags):
        if f:b[k]^=1;b[k+1]^=1
    if shift:bp=rotate(bp,1)
    return bp,b,p,dict(waits=waits,flags=flags,terminal=terminal,shift=shift,stash=stash)


def truth(s):
    n=len(s);P=[0]*n+[1];R=trim(s);bp=[];b=[1]
    while R:
        q,r=divide(P,R);nb=plus(bp,product(q,b));P,R=R,r;bp,b=b,nb
        if len(R)<=n//2 or not R:break
    return monic(b)


def evaluate(poly,x):
    r=0
    for v in reversed(poly):r=mul(r,x)^v
    return r


if __name__=='__main__':
    import json
    from pathlib import Path
    cases=list(map(json.loads,Path('parameters/exception_cases_1m.jsonl').read_text().splitlines()))
    groups={};wrong=[]
    for i,c in enumerate(cases):
        bp,b,p,meta=decode(c['syndromes']);bad=trim(b)!=truth(c['syndromes'])
        for pos,val in zip(c['positions'],c['values']):
            x=EX[-pos%4095];der=sum([]) if False else 0
            for j in range(1,len(b),2):der^=mul(b[j],EX[(-pos*(j-1))%4095])
            got=mul(EX[-pos*70%4095],inv(mul(mul(p,evaluate(bp,x)),der)))
            bad|=got!=val
        if bad:wrong.append(dict(index=i,shot=c['shot'],stage=c['first_stage'],gap_count=c['gap_count'],meta=meta))
        groups[c['first_stage']]=groups.get(c['first_stage'],0)+1
    report=dict(status='PASS' if not wrong else 'RESIDUAL_FAILURES',cases=len(cases),wrong=len(wrong),failures=wrong,stages=groups)
    Path('results/serial_reference_exceptions.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='failures'}),flush=True)
