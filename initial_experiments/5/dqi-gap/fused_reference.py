"""Compare proposed blocked local Schur coordinates with full serial EEA state."""
import json
from pathlib import Path
from classical_repair import mul,inv
from check_gap_block import sp,sd
from check_gap_serial import reference


def moments(bp,b,c,p,stash,carry,shift,count,k):
    raw=list(reversed(c))[:2*count]
    if carry:
        f=stash+raw;n=len(f)
        d=[] if k==1 else bp[::-1]
        den=[inv(p),0]+sp(f,d,max(0,n-2))
        h=sd(sp(f,b[::-1],n),den,n)
        return h[:2],h[2:]
    n=len(raw);d=bp[::-1];skip=3 if shift else 2
    den=[inv(p)]+[0]*(skip-1)+sp(raw,d,max(0,n-skip))
    return list(stash),sd(sp(raw,b[::-1],n),den,n)


def block_check(state,k,count,waits,flags,terminal):
    bp,b,c,p,stash=state;carry=waits[k-1];shift=(waits[k-2] if k>=2 else 0)^terminal
    hs,h=moments(bp,b,c,p,stash,carry,shift,count,k)
    for j in range(count):
        kk=k+j;carry=waits[-1];shift=(waits[kk-2] if kk>=2 else 0)^terminal
        view=h[2*j:]
        if carry:local=(hs if j==0 else h[2*j-2:2*j])+view
        else:local=view
        actuals,actual=moments(bp,b,c,p,stash,carry,shift,1,kk)
        actualfull=actuals+actual if carry else actual
        assert local[:len(actualfull)]==actualfull,dict(k=kk,block=k,j=j,carry=carry,shift=shift,local=local[:len(actualfull)],actual=actualfull)
        wait=int(not(view[0] or not view[1] or carry or shift))
        complete=0
        if not wait and not carry:
            complete=int(not view[0]);view[0]^=complete
            beta=view[0];out=sd([beta],view,len(view));view[1:]=out[1:]
        elif carry:
            d=local[1];u=sd([1],[1]+[mul(x,inv(d)) for x in local[2:]],len(local)-1)
            local[2:]=u[1:];view=local[2:]
        h[2*j:]=view
        bp,b,c,p,stash,flag=reference(bp,b,c,p,stash,carry,wait,shift)
        assert flag==complete
        terminal^=(waits[kk-2] if kk>=2 else 0)&flag;waits.append(wait);flags.append(flag)
    return (bp,b,c,p,stash),terminal


if __name__=='__main__':
    cases=list(map(json.loads,Path('parameters/exception_cases_1m.jsonl').read_text().splitlines()))[:128]
    for ci,case in enumerate(cases):
        s=case['syndromes'];wait=int(not s[-1] and s[-2]);waits=[wait];flags=[int(not s[-1] and not wait)];terminal=0
        if wait:state=([1],[0,1],s[:-2],1,[0,s[-2]])
        else:
            s=s.copy();s[-1]^=flags[0];p=inv(s[-1]);al=mul(s[-2],p)
            state=([1],[al,1],[x^mul(al,y) for x,y in zip(s[:-2],s[1:-1])],p,[0,0])
        for k in range(1,35,3):state,terminal=block_check(state,k,min(3,35-k),waits,flags,terminal)
    print(json.dumps(dict(status='PASS',cases=len(cases))),flush=True)
