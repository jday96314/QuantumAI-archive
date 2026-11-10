"""Model the syndrome-consuming BM template, including inactive borrowed slots."""
import argparse,json,random
from probe_wronskian import p,mul,evaluate
inv=p.ginv

def model(original,offset=3):
    S=list(original);sig=[1]+[0]*35;extra=[1]+[0]*(offset-1);limit={3:67,4:65,5:63}[offset]
    park=0;parkindex=0;filled=False;active=True;d=0;stack=False;reason=None;materialized=set()
    def get(k):return extra[k] if k<offset else S[k-offset]
    def put(k,value):
        if k<offset:extra[k]=value
        else:S[k-offset]=value
    def stop(why):
        nonlocal active,reason
        if active:reason=why;active=False
    def rotate(size):
        for k in reversed(range(1,size)):
            a,c=get(k),get(k-1);put(k,c);put(k-1,a)
    for r in range(70):
        size=min(37,max(p.KMSV[r],p.KMMV[r])+2)
        deferred=offset==3 and r%2 and r<=limit
        if not deferred:materialized.update(range(size));rotate(size)
        if r and r<69:
            guard=size-1 if deferred else 0
            if guard in materialized and get(guard) and active:stop(('bank',r))
            if d<0 and stack and active:stop(('stack',r))
        if r%2 and r<=limit:
            l=(r+1)//2;t=l-1
            if active and filled and sig[l]==0:stop(('secondpark',r))
            if active and sig[l]==0:
                S[t],park=park,S[t];filled=not filled;parkindex=l
            S[t]=mul(S[t],sig[l] or 1)
            if active:
                for i in range(p.IMAXV[r]+1):
                    if i!=l:S[t]^=mul(sig[i],S[r-i])
                D=S[t];S[t]=0
            else:D=0
            if parkindex==l and active and D==0:stop(('done',r))
        else:
            restore=r%2==0 and r>0 and parkindex==r//2 and active
            t=r//2-1
            if restore:S[t],park=park,S[t]
            D=0
            for i in range(p.IMAXV[r]+1):D^=mul(sig[i],S[r-i])
            if restore:S[t],park=park,S[t]
        if deferred:materialized.update(range(size));rotate(size)
        ctl=active and D!=0 and d>=0 and r!=69;q=mul(D,sig[0]);m=(inv(q) if ctl else q) if active else 0
        if ctl:
            for k in range(p.KMSV[r]+1):
                if r==1 and k==1:continue
                old=sig[k];sig[k]=get(k);put(k,old)
        for k in range(1 if r==69 else 0,p.KMMV[r]+1):sig[k]^=mul(m,get(k))
        if active and d<0:stack=True
        # Counter proxy: one larger than actual dW after round68, for a simpler degree formula.
        if r<69:d=-d if ctl else d+1
    return sig,[get(k) for k in range(37)],d,S,reason

def normalize_bank(bank,L):
    bank=list(bank)
    if L<=34:
        # Completion freezes the algebra, but the prefix rotations continue.
        # Stop just before the completion-round rotation (r=2L+1).
        for r in range(69,2*L,-1):
            size=min(37,max(p.KMSV[r],p.KMMV[r])+2)
            bank[:size]=bank[1:size]+bank[:1]
    else:
        bank=bank[-1:]+bank[:-1]
    return bank

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--shots',type=int,default=1000);a=parser.parse_args()
    rng=random.Random(51927);stats={};failures=[]
    for shot in range(a.shots):
        L=shot%36;pos=rng.sample(range(4095),L);vals=[rng.randrange(1,4096) for _ in pos];S=[]
        for k in range(1,71):
            s=0
            for j,v in zip(pos,vals):s^=mul(v,p.EXP[k*j%4095])
            S.append(s)
        sig,bank,d,syn,reason=model(S);bank=normalize_bank(bank,L)
        siggood=all(evaluate(sig,p.EXP[-j%4095])==0 for j in pos)
        # Only the degree-bounded template participates in the proposed map.
        trimmed=bank[:L+2]+[0]*(35-L);H=[0]*36
        for i,x in enumerate(sig):
            for j,y in enumerate(trimmed):
                if (i+j)&1:H[(i+j-1)//2]^=mul(x,y)
        H=[p.EXP[p.LOG[x]*2048%4095] if x else 0 for x in H]
        good=siggood and bank[0]==0
        for j,v in zip(pos,vals):
            h=mul(evaluate(H,p.EXP[-j%4095]),p.EXP[L*j%4095]);good &= mul(v,mul(h,h))==1
        key=str(reason[0] if reason else None);stats[key]=stats.get(key,0)+1
        if not good:failures.append(dict(shot=shot,L=L,reason=reason,sigma_correct=siggood,constant=bank[0]))
    print(json.dumps(dict(shots=a.shots,scope='classical stream model, degree-trimmed template',stops=stats,failures=failures),indent=2))
