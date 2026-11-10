"""Independent exact decoding, tangent, Wronskian and C-fraction checks."""
from pathlib import Path
import hashlib,json,random
from algebra import Field,trim,sylvester_reference

HERE=Path(__file__).resolve().parent


def solve(F,A,b):
    a=[list(row)+[x] for row,x in zip(A,b)];n=len(b)
    for j in range(n):
        k=next(k for k in range(j,n) if a[k][j]);a[j],a[k]=a[k],a[j]
        a[j]=[F.mul(v,F.inv(a[j][j])) for v in a[j]]
        for i in range(n):
            if i!=j:
                q=a[i][j]
                a[i]=[x^F.mul(q,y) for x,y in zip(a[i],a[j])]
    return [row[-1] for row in a]


def tangent_locator(F,s,d,b):
    ds=[F.mul(s[i//2],s[i//2]) if not i%2 else 0 for i in range(len(s))]
    A=[[s[k-j] for j in range(d)] for k in range(d,2*d)]
    rhs=[]
    for k in range(d,2*d):
        v=ds[k-d]
        for j in range(d):v^=F.mul(ds[k-j],b[j])
        rhs.append(v)
    return trim(solve(F,A,rhs))


def peel(F,f):
    coefficients=[];f=list(f)
    while f and any(f):
        a=f[0]
        if not a:return None
        coefficients.append(a)
        invf=F.reciprocal(f,len(f))
        f=[F.mul(a,x) for x in invf[1:]]
    coefficients += [0]*len(f)
    return coefficients


def synthesize(F,coefficients):
    p,q=[],[1]
    for a in reversed(coefficients):p,q=F.scale(q,a),F.add(q,[0]+p)
    return p,q


def wronskian_root(F,a,b):
    W=F.add(F.product(a,F.derivative(b)),F.product(F.derivative(a),b))
    assert all(not v for v in W[1::2])
    return trim([F.sqrt(x) for x in W[::2]])


def main():
    salt='core-phase-algebra-independent-20261024-a'
    rng=random.Random(int.from_bytes(hashlib.sha256(salt.encode()).digest(),'little'))
    cases=roots=mask_checks=tangents=cf_cases=cf_gaps=0;degree_patterns={};gap_cases=0
    for F,maxd,repeats in [(Field(4,0x13),7,30),(Field(),35,4)]:
        for d in range(maxd+1):
            for trial in range(repeats):
                N=2*maxd
                positions=rng.sample(range(1,F.size),d)
                values=[rng.randrange(1,F.size) for _ in positions]
                # Force a missing top syndrome using only this independent
                # error constructor. Decoders still receive syndromes alone.
                if d>=2 and trial==0:
                    z=0
                    for a,e in zip(positions[:-1],values[:-1]):z^=F.mul(e,F.pow(a,N))
                    last=F.mul(z,F.inv(F.pow(positions[-1],N)))
                    if last:values[-1]=last
                s=F.syndrome(positions,values,N)
                b,omega,pattern=F.pade(s)
                assert b==F.locator(positions)==F.bm(s),(F.degree,d,trial,'locator')
                assert trim(F.product(b,s)[:N])==omega
                degree_patterns[str(pattern)]=degree_patterns.get(str(pattern),0)+1
                gap_cases+=int(any(x>1 for x in pattern))
                if d:
                    db=tangent_locator(F,s,d,b)
                    assert db==omega,(F.degree,d,'Riccati tangent')
                    tangents+=1
                cf=peel(F,s)
                if cf is None:cf_gaps+=1
                else:
                    pp,qq=synthesize(F,cf)
                    assert F.product(pp,b)==F.product(omega,qq)
                    cf_cases+=1
                # Construct A from exact compact identity only for checking
                # the new representation, never to produce a decoded answer.
                if d:
                    # Extended Euclid: inverse omega modulo B, then x^N.
                    r0,r1=b,omega;u0,u1=[],[1]
                    while r1:
                        quotient,r=F.divide(r0,r1)
                        r0,r1,u0,u1=r1,r,u1,F.add(u0,F.product(quotient,u1))
                    assert len(r0)==1
                    oi=F.scale(u0,F.inv(r0[0]))
                    A=F.divide(F.product([0]*N+[1],oi),b)[1]
                    H=wronskian_root(F,A,b)
                    # Direct even/odd Sylvester form, including arbitrary
                    # degree gaps in E/O; direct multiplication is the oracle.
                    C=[F.sqrt(x) for x in A[::2]];D=[F.sqrt(x) for x in A[1::2]]
                    E=[F.sqrt(x) for x in b[::2]];O=[F.sqrt(x) for x in b[1::2]]
                    assert H==F.add(F.product(C,O),F.product(D,E))
                    for a,e in zip(positions,values):
                        beta=F.inv(a);hv=F.evaluate(H,beta)
                        got=F.mul(F.pow(beta,N),F.inv(F.mul(hv,hv)))
                        assert got==e,(F.degree,d,'H phase value')
                        for _ in range(3):
                            mask=rng.randrange(F.size)
                            # The runtime coordinate mask is transformed by
                            # x -> beta^N*x^2, an invertible binary linear map.
                            transformed=sum(((mask & F.mul(F.pow(beta,N),F.mul(1<<j,1<<j))).bit_count()%2)<<j for j in range(F.degree))
                            assert (mask&e).bit_count()%2==(transformed&F.inv(hv)).bit_count()%2
                            mask_checks+=1
                        roots+=1
                cases+=1
    # General Wronskian update: no restriction to linear quotients.
    recurrences=0;sylvester_cases=0
    F=Field(4,0x13)
    for _ in range(300):
        A=[rng.randrange(16) for _ in range(5)];B=[rng.randrange(16) for _ in range(6)]+[1]
        q=[rng.randrange(16) for _ in range(rng.randrange(1,6))];lam=rng.randrange(1,16)
        new=F.add(F.product(q,B),F.scale(A,lam))
        lhs=wronskian_root(F,B,new)
        rhs=F.add(F.product([F.sqrt(x) for x in q[1::2]],B),F.scale(wronskian_root(F,A,B),F.sqrt(lam)))
        assert lhs==rhs;recurrences+=1
        n=rng.randrange(1,6);m=rng.randrange(n,n+3)
        f=[rng.randrange(16) for _ in range(m)]+[rng.randrange(1,16)]
        g=[rng.randrange(16) for _ in range(n)]+[rng.randrange(1,16)]
        c=[rng.randrange(16) for _ in range(n)];dd=[rng.randrange(16) for _ in range(m)]
        try:h=sylvester_reference(F,f,g,c,dd)
        except ZeroDivisionError:continue
        want=F.add(F.product(c,f),F.product(dd,g));want += [0]*(m+n-len(want))
        assert h==want;sylvester_cases+=1
    # The Riccati tangent is localized in the first TWO C-fraction heads.
    cf_tangents=0;attempts=0
    while cf_tangents<40:
        attempts+=1;f=[rng.randrange(1,16)]+[rng.randrange(16) for _ in range(15)]
        df=F.product(f,f)[:len(f)];df += [0]*(len(f)-len(df));heads=[];dheads=[]
        valid=True
        while f:
            if not f[0]:valid=False;break
            a,da=f[0],df[0];iv=F.reciprocal(f,len(f))
            iv2=F.product(iv,iv)[:len(f)];div=F.product(df,iv2)[:len(f)]
            div += [0]*(len(f)-len(div));dg=F.add(F.scale(iv,da),F.scale(div,a));dg += [0]*(len(f)-len(dg))
            assert dg[0]==0
            heads.append(a);dheads.append(da)
            f=[F.mul(a,x) for x in iv[1:]];df=dg[1:]
        if not valid:continue
        expected=[F.mul(heads[0],heads[0]),F.mul(heads[0],heads[1])]+[0]*(len(heads)-2)
        assert dheads==expected;cf_tangents+=1
    result=dict(status='PASS',salt=salt,rs_cases=cases,root_value_checks=roots,runtime_mask_checks=mask_checks,
                exact_tangent_checks=tangents,rs_cases_with_degree_gaps=gap_cases,degree_patterns=degree_patterns,
                c_fraction_reconstructed_cases=cf_cases,c_fraction_zero_head_charts_unimplemented=cf_gaps,
                general_quotient_wronskian_recurrences=recurrences,arbitrary_gap_sylvester_cases=sylvester_cases,
                localized_C_fraction_tangent_checks=cf_tangents,localized_tangent_attempts=attempts,
                scope='independent exact finite-field EEA and BM vs constructed errors; no native task-cost claim')
    (HERE/'results/algebra.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
