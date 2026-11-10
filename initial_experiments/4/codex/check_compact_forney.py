"""Check locator and inverse-square magnitude from compact Euclid state.

All weights are exercised. Degree-gap cases are counted and skipped; the
diagnostic is not a native complete-decoder verification.
"""
import argparse,json,random,time
from pathlib import Path
from compact_euclid import mul,inv,product,xor
from monic_compact_euclid import first,step


def evaluate(poly,x):
    value=0
    for coefficient in reversed(poly):value=mul(value,x)^coefficient
    return value


def derivative(poly):return [poly[j] if j%2 else 0 for j in range(1,len(poly))]


def test(N,trials):
    from import_published import load
    m=load();rng=random.Random(104713);begin=time.time();checks=0;gaps=0;weights={};roots=0
    for weight in range(1,N//2+1):
        for _ in range(trials):
            positions=rng.sample(range(1,4096),weight);values=[rng.randrange(1,4096) for x in positions]
            S=[0]*N
            for x,e in zip(positions,values):
                power=e
                for j in range(N):S[j]^=power;power=mul(power,x)
            try:
                state=first(S)
                for k in range(1,weight):state=step(k,*state)
            except ValueError:gaps+=1;continue
            bp,b,c,pinv=state
            assert not any(c),(weight,'completion')
            expected=[1]
            for position in positions:expected=product(expected,[inv(position),1])
            assert b==expected,(weight,'locator')
            wronskian=xor(product(derivative(b),bp),product(b,derivative(bp)))
            assert not any(wronskian[1::2]),(weight,'odd Wronskian')
            # Folding sqrt(pinv) into H makes the magnitude only an inverse
            # square times the known field position power x^(N-1).
            H=[m.gpow(mul(wronskian[j],pinv),2048) for j in range(0,len(wronskian),2)]
            for position,e in zip(positions,values):
                root=inv(position);h=evaluate(H,root)
                predicted=mul(m.gpow(root,N-1),m.gpow(inv(h),2))
                assert predicted==e,(weight,'magnitude',predicted,e)
                assert evaluate(b,root)==0;roots+=1
            checks+=1;weights[str(weight)]=weights.get(str(weight),0)+1
    out=dict(status='PASS',scope='classical compact-state exact locator, completion and inverse-square magnitude, all nonzero weights, skipped generic degree gaps',N=N,trials_per_weight=trials,cases=checks,error_roots=roots,skipped_degree_gaps=gaps,weights=weights,seconds=time.time()-begin)
    Path(f'results/compact_forney_n{N}_check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=70);p.add_argument('--trials',type=int,default=20);a=p.parse_args();test(a.n,a.trials)
