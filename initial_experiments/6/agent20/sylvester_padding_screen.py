"""Independent promised-input screen for the existing armed Sylvester map."""
import json
from pathlib import Path
import numpy as np
from import_published import load

m=load()
def prod(a,b,n):
    c=[0]*n
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<n:c[i+j]^=m.gmul(x,y)
    return c

def sylvester(f,g,c,d,repair_state=None):
    f=list(f);g=list(g);c=list(c);d=list(d)
    if not g[-1]:g[-1]=1
    n=len(g)-1;mm=len(f)-1
    if not n:return [m.gmul(v,g[0]) for v in d]
    for i in reversed(range(n,mm+1)):
        f[i]=m.gmul(f[i],m.ginv(g[-1]))
        for j in range(n):f[i-n+j]^=m.gmul(f[i],g[j])
    q,r=f[n:],f[:n];p=prod(q,c,mm);d=[a^b for a,b in zip(d,p)]
    k=n-1;actual_r=list(r)
    repair=bool(repair_state is not None and not repair_state[0] and len(r)>1 and not r[-1] and r[-2])
    gamma=0
    padding_floor=repair_state[1] if repair_state is not None and len(repair_state)>1 else 0
    if repair:
        repair_state[0]=True
        gamma=c[padding_floor]
        c=c[:padding_floor]+c[padding_floor+1:]+[0]
        r=[0,*r[:-1]]
    low=sylvester(g,r,d[:k],c,repair_state);high=d[k:]
    if repair:
        for i,x in enumerate(actual_r[:-1]):low[i+padding_floor]^=m.gmul(gamma,x)
    cross=prod(g,high,len(g)+len(high)-1)
    for i in range(n):low[i+k]^=cross[i]
    for j in range(len(high)):
        high[j]=m.gmul(high[j],g[-1])
        for ell in range(j+1,len(high)):
            idx=n+j-ell
            if idx>=0:high[j]^=m.gmul(g[idx],high[ell])
    return low+high

def run(n=35,shots=100,repair=False):
    assert n%2;t=n//2;rng=np.random.default_rng(2026110401);rows=[]
    for degree in range(n+1):
        good=0;bad=0
        for _ in range(shots):
            b=[0]*(n-degree)+rng.integers(0,4096,degree).tolist()+[1]
            a=[0]*(n-degree)+rng.integers(0,4096,degree).tolist()
            a+=[0]*(n-len(a))
            sqrt=lambda x:m.gpow(x,2048)
            E=list(map(sqrt,b[::2]));O=list(map(sqrt,b[1::2]));C=list(map(sqrt,a[::2]));D=list(map(sqrt,a[1::2]))
            h=sylvester(E,O,D,C[:-1],[False,(n-degree)//2] if repair else None)+[C[-1]]
            for i,x in enumerate(O[:-1]):h[t+i]^=m.gmul(C[-1],x)
            expected=[x^y for x,y in zip(prod(C,O,n),prod(D,E,n))]
            if h==expected:good+=1
            else:bad+=1
        row=dict(degree=degree,padding=n-degree,good=good,bad=bad);rows.append(row);print(json.dumps(row),flush=True)
    Path('results/sylvester_padding'+('_repair_v2' if repair else '_screen')+'.json').write_text(json.dumps(rows,indent=2)+'\n')

if __name__=='__main__':run()
