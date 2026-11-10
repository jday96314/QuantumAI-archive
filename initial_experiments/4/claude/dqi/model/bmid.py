"""Classical check: at BM step k, is the reversed-window evaluator inverse to the template mod sigma*?"""
import random
import kickmix as km
F=km.GF2Field(12);mul=F.mul;inv=F.invert;alpha=F.primitive_element
def pmul(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):r[i+j]^=mul(x,y)
    return r
def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return a
def pmod(a,m):
    a=trim(a);m=trim(m);il=inv(m[-1])
    while len(a)>=len(m) and any(a):
        c=mul(a[-1],il);s=len(a)-len(m)
        for i,y in enumerate(m):a[s+i]^=mul(c,y)
        a=trim(a)
    return a
rng=random.Random(1)
v=30;pos=rng.sample(range(4095),v);val=[rng.randrange(1,4096) for _ in pos]
S=[0]*70
for p,e in zip(pos,val):
    for i in range(70):S[i]^=mul(e,F.pow(alpha,((i+1)*p)%4095))
# standard BM
C=[1];Bp=[1];L=0;m=1;b=1
for n in range(70):
    d=S[n]
    for i in range(1,L+1):
        if i<len(C):d^=mul(C[i],S[n-i])
    if d==0:m+=1;continue
    coef=mul(d,inv(b));T=list(C)
    shifted=[0]*m+[mul(coef,x) for x in Bp]
    C=[ (C[i] if i<len(C) else 0)^(shifted[i] if i<len(shifted) else 0) for i in range(max(len(C),len(shifted)))]
    if 2*L<=n:L=n+1-L;Bp=T;b=d;m=1
    else:m+=1
    k=n+1
    if k in (20,30,40,47,60,70):
        sig=trim(C);tau=[0]*m+Bp  # template aligned (x^m B)
        Ls=len(sig)-1
        W=S[k-Ls:k]  # window
        # reversed polys
        om0=pmul(sig,S[:k])[:Ls]
        th0=pmul(tau,S[:k])[:k]
        lhs=[x^y for x,y in zip(pmul(sig,th0)+[0]*200, pmul(tau,om0)+[0]*200)]
        nz=[i for i,x in enumerate(trim(lhs)) if x]
        print('berlekamp nonzero degrees',nz[:5],nz[-3:] if nz else None)
        sigr=list(reversed(sig+[0]*(Ls+1-len(sig))))
        Wr=list(reversed(W))
        om=pmul(sigr,Wr)[:Ls]
        taur=list(reversed(tau))
        prod=pmod(pmul(om,taur),sigr)
        print(k,L,len(tau)-1,'prod mod sigma*:',prod[:4],'len',len(trim(prod)))
