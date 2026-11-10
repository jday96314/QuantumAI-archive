"""Check actual round-4 cofactor conventions against independent errors.

Calls the read-only round-4 classical state model; no recurrence is copied,
implemented or modified. Generic degree gaps are reported, not repaired.
The independent oracle constructs syndromes directly from known errors.
"""
import json
import random
from pathlib import Path
from import_published import load
from published_dqi import isomorphism
from monic_compact_euclid import first,step
from oracle import known_error_polynomials,evaluate,derivative,mul,inv,power


def main():
    rng=random.Random(546117);m=load();fw,bk,transform=isomorphism(m)
    cases=roots=0;gaps=[];all_degrees={}
    for L in range(1,36):
        good=0
        for trial in range(20):
            pos=rng.sample(range(4095),L);values=[rng.randrange(1,4096) for _ in pos]
            sigma,omega,S=known_error_polynomials(pos,values)
            S_pub=[transform(v,fw) for v in S]
            try:
                state=first(S_pub)
                for k in range(1,L):state=step(k,*state)
            except ValueError as e:
                gaps.append(dict(degree=L,trial=trial,reason=str(e)));continue
            bp,b,c,pinv=state
            assert not any(c),(L,'completion')
            B=[transform(v,bk) for v in b]
            A=[transform(m.gmul(v,pinv),bk) for v in bp]
            derivativeB=derivative(B)
            monic=[mul(v,inv(sigma[-1])) for v in sigma]
            assert B==monic,(L,'locator orientation')
            for j,e in zip(pos,values):
                beta=power(2,-j)
                assert evaluate(B,beta)==0
                got=mul(power(beta,70),inv(mul(evaluate(A,beta),evaluate(derivativeB,beta))))
                assert got==e,(L,j,got,e)
                roots+=1
            good+=1;cases+=1
        all_degrees[L]=good
    out=dict(status='PASS',scope='read-only round-4 cofactor model versus independent known-error oracle; generic paths only, not native decoder acceptance',
             seed=546117,cases=cases,root_checks=roots,generic_degree_gaps=gaps,cases_per_degree=all_degrees,
             polynomial_layout='B monic, A=pinv*Bprev; beta=alpha^(-j); e=beta^70/(A(beta)*Bprime(beta))',
             evaluator_conversion='omega_monic = x^70 * inverse(A) mod B; F=x^70*inverse(A)*inverse(Bprime) mod B')
    Path('results/control_identity.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':main()
