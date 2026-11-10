"""Classical packed-array model of one repaired degree gap in fixed EEA."""
import numpy as np
import kickmix as km

GF = km.GF2Field(12)


def mul(a,b):
    return GF.mul(int(a),int(b))


def div(a,b):
    return mul(a, GF.pow(int(b),4094))


def decode(syndrome, trace=False):
    T = len(syndrome)
    t = T//2
    A = [0]*(T+1)
    B = list(map(int,syndrome))+[1]
    A[T] = 1
    counter = 0
    index = 0
    for i in range(1,t+1):
        X,Y = (A,B) if i%2 else (B,A)
        active = counter == i-1 and Y[T-i] != 0
        bad = counter == i-1 and Y[T-i] == 0 and Y[T-i-1] != 0 and index == 0
        pending = index == i-1 and index != 0
        if active:
            q2 = div(X[T+2-i],Y[T-i]) if pending else 0
            q1 = div(X[T+1-i] ^ mul(q2,Y[T-i-1]),Y[T-i])
            q0 = div(X[T-i] ^ mul(q1,Y[T-i-1]) ^ mul(q2,Y[T-i-2]),Y[T-i])
            qs = [q0,q1,q2]
            for k in range(min(T+1,T-i+3)):
                for d,q in enumerate(qs):
                    if 0 <= k-d <= T-i:
                        X[k] ^= mul(q,Y[k-d])
            for k in range(i+2):
                for d,q in enumerate(qs):
                    if 0 <= k-d <= i-1:
                        X[T-k] ^= mul(q,Y[T-(k-d)])
            ys = [Y[T-i+1+int(pending)+j] if T-i+1+int(pending)+j <= T else 0 for j in range(3)]
            h = div(X[T-i],ys[0])
            m = div(X[T-i+1] ^ mul(h,ys[1]),ys[0])
            lo = div(X[T-i+2] ^ mul(m,ys[1]) ^ mul(h,ys[2]),ys[0]) if pending else 0
            recovered = [lo,m,h] if pending else [m,h,0]
            assert recovered == qs, (i,index,qs,recovered,ys)
        if bad:
            A,B = B,A
            index = i
        counter += int(active or bad)
        if trace:
            print(i,'active',active,'skip',bad,'counter',counter,'index',index)
    if counter%2:
        A,B = B,A
    return B,index,counter


def evaluate(coeffs,x):
    out = 0
    for c in reversed(coeffs):
        out = mul(out,x)^int(c)
    return out


def decode_offset(syndrome, trace=False):
    T=len(syndrome);t=T//2
    A=[0]*(T+1);B=list(map(int,syndrome))+[1];A[T]=1
    counter=0;index=0
    for i in range(1,t+1):
        X,Y=(A,B) if i%2 else (B,A)
        off=int(index!=0)
        bad=i<t and counter==i-1 and not off and Y[T-i]==0 and Y[T-i-1]!=0
        active=i+off<=t and counter==i-1 and (Y[T-i-off]!=0 or bad)
        if active:
            ytop=T-i-off-int(bad);xtop=T+1-i-off
            h=div(X[xtop],Y[ytop])
            m=div(X[xtop-1]^mul(h,Y[ytop-1]),Y[ytop])
            lo=div(X[xtop-2]^mul(m,Y[ytop-1])^mul(h,Y[ytop-2]),Y[ytop]) if bad else 0
            qs=[lo,m,h] if bad else [m,h,0]
            rY=Y[:T-i-off+1]
            sY=[Y[T-d] for d in range(i+off)]
            for k in range(T+1):
                for d,q in enumerate(qs):
                    if 0<=k-d<len(rY):X[k]^=mul(q,rY[k-d])
            for k in range(i+off+int(bad)+1):
                for d,q in enumerate(qs):
                    if 0<=k-d<len(sY):X[T-k]^=mul(q,sY[k-d])
            xstop=T-i-off-int(bad);ystop=T-i+1-off
            h=div(X[xstop],Y[ystop])
            yn=Y[ystop+1] if ystop+1<=T else 0
            y2=Y[ystop+2] if ystop+2<=T else 0
            m=div(X[xstop+1]^mul(h,yn),Y[ystop])
            lo=div(X[xstop+2]^mul(m,yn)^mul(h,y2),Y[ystop]) if bad else 0
            got=[lo,m,h] if bad else [m,h,0]
            assert got==qs,(i,index,bad,got,qs)
        if bad:index=i
        counter+=int(active)
    if counter%2:A,B=B,A
    return B,index,counter


def test(decoder=decode):
    rng = np.random.default_rng(32)
    for T in [6,10,20,70]:
        for forced in [False,True]:
            for shot in range(100):
                w = int(rng.integers(2,T//2+1))
                pos = rng.choice(4095,w,replace=False)
                roots = [GF.pow(GF.primitive_element,int(p)) for p in pos]
                values = list(map(int,rng.integers(1,4096,w)))
                if forced:
                    top = 0
                    for e,r in zip(values[:-1],roots[:-1]):
                        top ^= mul(e,GF.pow(r,T))
                    values[-1] = div(top,GF.pow(roots[-1],T))
                synd = [0]*T
                for e,r in zip(values,roots):
                    for j in range(T):
                        synd[j] ^= mul(e,GF.pow(r,j+1))
                result,index,count = decoder(synd)
                sig = list(reversed(result[T-T//2:T+1]))
                om = result[:T//2]
                for p,r,e in zip(pos,roots,values):
                    x = GF.invert(r)
                    ss = evaluate(sig,x)
                    dd = evaluate([sig[j] if j%2 else 0 for j in range(len(sig))],x)
                    ee = div(mul(x,evaluate(om,x)),dd)
                    assert ss == 0 and ee == e,(T,forced,shot,p,e,ee,ss,index,count)
            print('PASS',T,'forced first degree gap',forced,flush=True)


if __name__ == '__main__':
    import sys
    test(decode_offset if len(sys.argv)>1 else decode)
