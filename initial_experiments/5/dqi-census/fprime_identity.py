"""Independent error-polynomial check for A -> (A B')^-1 and mask folding."""
import json
from pathlib import Path
import numpy as np
from oracle import product, divide, inverse_mod, derivative, evaluate, mul, power, MASK_TO_DUAL, trace
ROOT=Path(__file__).resolve().parent

def main():
    runs=[]
    for seed in [2026102427,2026102428]:
        rng=np.random.default_rng(seed);cases=roots_checked=masks_checked=0
        for weight in range(36):
            for trial in range(4):
                positions=rng.choice(4095,weight,replace=False)
                beta=[power(2,(-int(j))%4095) for j in positions]
                e=rng.integers(1,4096,weight).tolist();B=[1]
                for x in beta:B=product(B,[x,1])
                # Omega=sum e*B/(x+beta), independently from geometric syndromes.
                omega=[0]*max(1,weight)
                for x,value in zip(beta,e):
                    q,r=divide(B,[x,1]);assert not r
                    for i,c in enumerate(q):omega[i]^=mul(value,c)
                if weight:
                    A=divide(product([0]*70+[1],inverse_mod(omega,B)),B)[1]
                    F= inverse_mod(divide(product(A,derivative(B)),B)[1],B)
                    back=inverse_mod(divide(product(F,derivative(B)),B)[1],B)
                    assert back==A
                    reconstructed=[mul(power(x,70),evaluate(F,x)) for x in beta]
                    assert reconstructed==e
                else:F=[]
                masks=rng.integers(0,4096,4095).tolist()
                direct=folded=0
                for j,x,value in zip(positions,beta,e):
                    direct^=(int(masks[int(j)])&value).bit_count()&1
                    mu=int(MASK_TO_DUAL[int(masks[int(j)])])
                    folded^=trace(mul(mul(mu,power(x,70)),evaluate(F,x)))
                assert direct==folded
                cases+=1;roots_checked+=weight;masks_checked+=len(masks)
        runs.append(dict(seed=seed,cases=cases,root_values=roots_checked,independent_mask_words=masks_checked,failures=0))
    out=dict(status='PASS',runs=runs,scope='classical independent polynomial identity in host basis, all weights0..35, two fixed seeds; no decoding recurrence and no circuit cost/failure claim')
    (ROOT/'results/fprime_identity.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
