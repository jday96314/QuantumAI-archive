"""Classical completion padding preserving real roots and recoverable H.

Pad B to degree35 using distinct linear factors. Bprev mixes with B so the
pair remains coprime unless an artificial root overlaps a real one. H can be
un-padded by affine polynomial divisions after a fixed-degree Wronskian map.
This is an identity experiment, not a quantum circuit or submission.
"""
import argparse,json,random,time
from pathlib import Path
from compact_euclid import mul,inv,product,xor,divide
from monic_compact_euclid import first
from block_euclid_ratio import series_product,series_divide
from check_compact_forney import evaluate,derivative
from import_published import load
M=load();OMEGA=next(x for x in range(2,4096) if mul(x,x)^x==1)


def scaled(poly,c):return [mul(x,c) for x in poly]


def pad(bp,b,c):
    # q=omega in GF4, 1+q=omega^2; both output cofactors stay monic.
    P=[c,1]
    return scaled(xor(product(P,bp),scaled(b,OMEGA)),inv(1^OMEGA)),product(P,b)


def trace(value):
    result=0
    for _ in range(12):result^=value;value=mul(value,value)
    return result


def quadratic_constants():
    # Cheap constant linear maps, with trace one so x^2+x+b is irreducible.
    candidates=[b for b in range(1,4096) if trace(b)==1]
    candidates.sort(key=lambda b:sum(mul(b,1<<i).bit_count() for i in range(12)))
    return candidates[:18]


QUADRATICS=quadratic_constants()


def H_of(bp,b):
    w=xor(product(derivative(b),bp),product(b,derivative(bp)))
    assert not any(w[1::2])
    return [M.gpow(x,2048) for x in w[::2]]


def unpad_b_h(b,H,roots):
    for c in reversed(roots):
        P=c if isinstance(c,list) else [c,1]
        b,remainder=divide(b,P);assert not remainder
        if len(P)==3:numerator=xor(H,b)
        else:numerator=xor(scaled(H,M.gpow(1^OMEGA,2048)),scaled(b,M.gpow(OMEGA,2048)))
        H,remainder=divide(numerator,P);assert not remainder
    return b,H


def block_step(k,count,bp,b,c,pinv,quadratic=False):
    n=2*count;oldb=list(b);oldc=list(c);t=list(reversed(c))[:n]
    f=series_divide(series_product(t,list(reversed(b)),n),[inv(pinv),0]+series_product(t,list(reversed(bp)),max(0,n-2)),n)
    flags=[];roots=[];collisions=0;quotients=[]
    for j in range(count):
        zero=not f[0];flags.append(zero)
        if zero:f[0]=1
        beta=f[0];u=series_divide([beta],f,len(f));alpha=u[1];f=u[2:]
        if zero and alpha:raise ValueError('degree gap')
        quotients.append((alpha,beta))
    j=0
    while j<count:
        alpha,beta=quotients[j]
        if flags[j]:
            if quadratic and (k+j)%2==1 and j+1<count and flags[j+1]:
                P=[QUADRATICS[(k+j-1)//2],1,1]
                bp,b=xor(product(P,bp),b),product(P,b);roots.append(P);j+=2;continue
            root=0 if quadratic else k+j+1
            collisions+=not evaluate(b,root);bp,b=pad(bp,b,root);roots.append(root)
        else:bp,b=b,xor(product([alpha,1],b),scaled(bp,beta))
        pinv=mul(pinv,inv(beta));j+=1
    numerator=series_product(list(reversed(b)),list(reversed(oldc)),len(oldc))
    c=list(reversed(series_divide(numerator[2*count:],list(reversed(oldb)),len(oldc)-2*count)))
    return (bp,b,c,pinv),flags,roots,collisions


def test(N,width,trials,quadratic=False):
    rng=random.Random(917413);begin=time.time();cases=0;gaps=0;collision_cases=0;error_roots=0
    for weight in range(1,N//2+1):
        for _ in range(trials):
            positions=rng.sample(range(1,4096),weight);errors=[rng.randrange(1,4096) for p in positions];S=[0]*N
            for position,error in zip(positions,errors):
                power=error
                for j in range(N):S[j]^=power;power=mul(power,position)
            try:
                state=first(S);flags=[False];roots=[];collisions=0;k=1
                while k<N//2:
                    count=min(width,N//2-k);state,ff,rr,cc=block_step(k,count,*state,quadratic=quadratic)
                    flags.extend(ff);roots.extend(rr);collisions+=cc;k+=count
            except ValueError:gaps+=1;continue
            assert flags==[False]*weight+[True]*(N//2-weight),(weight,'completion flags',flags)
            bp,b,c,pinv=state;assert not c
            assert len(b)==N//2+1 and len(bp)==N//2
            H=H_of(bp,b);assert len(H)==N//2 and H[-1]==1
            b,H=unpad_b_h(b,H,roots)
            expected=[1]
            for position in positions:expected=product(expected,[inv(position),1])
            assert b==expected,(weight,'locator')
            for position,error in zip(positions,errors):
                root=inv(position);h=evaluate(H,root)
                predicted=mul(mul(inv(pinv),M.gpow(root,N-1)),M.gpow(inv(h),2))
                assert predicted==error,(weight,'magnitude');error_roots+=1
            collision_cases+=bool(collisions);cases+=1
    out=dict(status='PASS',scope='classical padded block Euclid, completion flags, fixed degree monic cofactors, Wronskian unpadding, locator and inverse-square magnitude; generic paths only',N=N,block=width,quadratic=quadratic,trials_per_weight=trials,cases=cases,error_roots=error_roots,skipped_degree_gaps=gaps,artificial_root_collision_cases=collision_cases,omega=OMEGA,quadratic_constants=QUADRATICS if quadratic else [],seconds=time.time()-begin)
    tag=f'coprime_padding_n{N}_b{width}'+('_quadratic' if quadratic else '')
    Path(f'results/{tag}_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=70);p.add_argument('--block',type=int,default=4);p.add_argument('--trials',type=int,default=20);p.add_argument('--quadratic',action='store_true');a=p.parse_args();test(a.n,a.block,a.trials,a.quadratic)
