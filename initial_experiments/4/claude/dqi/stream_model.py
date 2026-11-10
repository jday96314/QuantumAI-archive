"""Classical model for syndrome-consuming q-form Berlekamp--Massey.

Only phase approximation is intended: the eventual quantum implementation
must restore every syndrome, including branches that stop early.
"""
import random,json
from import_published import load

pub=load();mul=pub.gmul;inv=pub.ginv

def decode(original,offset=5,trace=False):
    history=[]
    S=list(original);sig=[0]*36;sig[0]=1
    extra=[0]*offset
    def get(k):return extra[k] if k<offset else S[k-offset]
    def put(k,v):
        if k<offset:extra[k]=v
        else:S[k-offset]=v
    extra[0]=1
    active=True;park=None;parkindex=None;d=0;stack=False;aborted=None
    def stop(reason,r):
        nonlocal active,aborted
        active=False;aborted=(reason,r)
    for r in range(70):
        size=min(37,max(pub.KMSV[r],pub.KMMV[r])+2)
        if active:
            newindex=size-1
            if newindex>=offset:
                assert newindex-offset<(r//2), (r,size,newindex-offset)
        for j in reversed(range(1,size)):
            v=get(j);put(j,get(j-1));put(j-1,v)
        if active and r<69 and get(0):stop('bank',r)
        if active and r<69 and d<0 and stack:stop('stack',r)
        target=(r-1)//2 if r%2 else None
        if active and r%2:
            l=(r+1)//2
            if any(sig[l+1:]):stop('high_degree',r)
            elif sig[l]==0 and park is not None:stop('second_gap_or_done',r)
        if active and r%2 and sig[(r+1)//2]==0:
            parkindex=target;park=S[target];S[target]=0
        restore=(active and r%2==0 and parkindex is not None and r==2*(parkindex+1))
        if restore:
            assert S[parkindex]==0
            S[parkindex],park=park,S[parkindex]
        D=0
        for i in range(pub.IMAXV[r]+1):D^=mul(sig[i],S[r-i])
        if restore:S[parkindex],park=park,S[parkindex]
        if active and r%2:
            # Nonzero leading coefficient scales the old syndrome into D.
            # A zero leading coefficient parks that input first instead.
            S[target]=D
            if sig[(r+1)//2]==0 and D==0:stop('done',r)
        ctl=active and D!=0 and d>=0 and r!=69
        q=mul(D,sig[0])
        m=(inv(q) if ctl else q) if active else 0
        if ctl:
            for k in range(pub.KMSV[r]+1):
                if r==1 and k==1:continue
                tmp=sig[k];sig[k]=get(k);put(k,tmp)
        for k in range(1 if r==69 else 0,pub.KMMV[r]+1):sig[k]^=mul(m,get(k))
        if active:
            assert sig[0]!=0,(r,'zero sigma0')
            if d<0:stack=True
        if active and r%2:S[target]=0
        elif r%2 and aborted==('done',r):S[target]=0
        if r<69:d=(-d) if ctl else d+1
        if trace:history.append((sig.copy(),d))
    return (sig,aborted,history) if trace else (sig,aborted)

def samples(count=1000):
    rng=random.Random(891)
    stats={};errors=0
    for shot in range(count):
        w=rng.randrange(1,36);pos=rng.sample(range(4095),w);val=[rng.randrange(1,4096) for _ in pos]
        S=[0]*70;expected=[0]*36;expected[0]=1
        for p,v in zip(pos,val):
            root=pub.EXP[p]
            for i in range(70):S[i]^=mul(v,pub.EXP[((i+1)*p)%4095])
            for i in reversed(range(1,36)):expected[i]^=mul(expected[i-1],root)
        sig,why=decode(S)
        mon=[mul(v,inv(sig[0])) for v in sig]
        bad=mon!=expected;errors+=bad
        key=str(why[0] if why else None);stats[key]=stats.get(key,0)+1
        if bad and errors<=5:print('error',shot,w,why,flush=True)
    print(json.dumps(dict(shots=count,locator_errors=errors,stops=stats)),flush=True)

if __name__=='__main__':
    import sys
    samples(int(sys.argv[1]) if len(sys.argv)>1 else 1000)
