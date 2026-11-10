"""Exact quotient-block extraction from a compact monic Euclidean state.

With P B+R Bprev=p*x^N, put T=reverse(R), A=reverse(B), D=reverse(Bprev).
The local quotient series is T*A/(p+z^2*T*D). Its first 2b coefficients
determine b generic degree-one Euclidean steps, independent of the full degree.
This is a classical identity test, not a completed reversible decoder.
"""
import json,random
from pathlib import Path
from compact_euclid import mul,inv,product,xor


def trim(a):
    a=list(a)
    while a and not a[-1]:a.pop()
    return a


def series_product(a,b,n):
    return (product(a,b)+[0]*n)[:n]


def series_divide(a,b,n):
    if not b or not b[0]:raise ValueError('zero series pivot')
    out=[0]*n;bi=inv(b[0])
    for i in range(n):
        value=a[i] if i<len(a) else 0
        for j in range(1,min(i,len(b)-1)+1):value^=mul(b[j],out[i-j])
        out[i]=mul(value,bi)
    return out


def normalized_step(P,R,Bprev,B):
    if not P or not R or len(P)!=len(R)+1:raise ValueError('degree gap or completion')
    beta=mul(R[-1],inv(P[-1]))
    alpha=mul(P[-2],inv(P[-1]))^mul(R[-2] if len(R)>1 else 0,inv(R[-1]))
    newR=trim(xor([mul(beta,v) for v in P],product([alpha,1],R)))
    newB=trim(xor([mul(beta,v) for v in Bprev],product([alpha,1],B)))
    return list(R),newR,list(B),newB,(alpha,beta)


def quotient_block(P,R,b):
    steps=[]
    for _ in range(b):
        P,R,_,_,q=normalized_step(P,R,[],[1]);steps.append(q)
    return steps


def local_series(R,Bprev,B,p,b):
    n=2*b;T=list(reversed(R))[:n];A=list(reversed(B));D=list(reversed(Bprev))
    numerator=series_product(T,A,n)
    denominator=[p,0]+series_product(T,D,max(0,n-2))
    return series_divide(numerator,denominator[:n],n)


def test():
    rng=random.Random(491291);trials=80;checks=0;rejected=0;stages={}
    for trial in range(trials):
        N=70;S=[rng.randrange(4096) for _ in range(N)]
        P=[0]*N+[1];R=S[:];Bp=[];B=[1]
        for k in range(N//2):
            if not R or len(P)!=len(R)+1:
                rejected+=1;break
            if k in [0,1,4,8,16,24,30,34]:
                for b in [1,2,4,8,16]:
                    if k+b>N//2:continue
                    f=local_series(R,Bp,B,P[-1],b)
                    direct=series_divide(list(reversed(R)),list(reversed(P)),2*b)
                    assert f==direct,(trial,k,b,'ratio')
                    # Quotients of x^(2b) / reverse(f) share the required prefix.
                    localP=[0]*(2*b)+[1];localR=list(reversed(f))
                    try:
                        expected=quotient_block(P,R,b)
                        got=quotient_block(localP,localR,b)
                    except ValueError:
                        rejected+=1;continue
                    assert got==expected,(trial,k,b,'quotient block',got,expected)
                    checks+=1;key=f'k{k}_b{b}';stages[key]=stages.get(key,0)+1
            try:P,R,Bp,B,_=normalized_step(P,R,Bp,B)
            except ValueError:rejected+=1;break
    out=dict(status='PASS',scope='classical exact local-series and normalized quotient blocks; generic degree-one paths only',
             N=70,trials=trials,block_checks=checks,skipped_degree_gaps=rejected,stages=stages)
    Path('results/block_euclid_ratio_check.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)


if __name__=='__main__':test()
