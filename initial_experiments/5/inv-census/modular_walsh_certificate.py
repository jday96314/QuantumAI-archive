"""Exact check that whole-turn angle lifts cannot sparsify single-bit Walsh.

A full-rank eleven-row, ten-column Walsh submatrix has a one-dimensional
integer left nullspace. Its primitive vector v requires v dot f to be an
integer for a phase function f to lie in that subtorus. Check every column
set. Modular rank is exact: all minors are smaller than the prime. Cofactor
reconstruction is verified with integer products and a primitive gcd.
"""
import itertools,math,json
from pathlib import Path
import numpy as np
P=(1<<31)-1
H=np.array([[1-2*((s&m).bit_count()%2) for m in range(16)] for s in range(11)],dtype=np.int64)

def rank_mod(mat):
    a=[[int(t)%P for t in row] for row in mat];rank=0
    for j in range(len(a[0])):
        piv=next((i for i in range(rank,len(a)) if a[i][j]),None)
        if piv is None:continue
        a[rank],a[piv]=a[piv],a[rank];iv=pow(a[rank][j],-1,P)
        a[rank][j:]=[t*iv%P for t in a[rank][j:]]
        for i in range(rank+1,len(a)):
            f=a[i][j]
            if f:a[i][j:]=[(aa-f*bb)%P for aa,bb in zip(a[i][j:],a[rank][j:])]
        rank+=1
    return rank

if __name__=='__main__':
    hist={};ranks={};maxv=-1;example=None
    for cols in itertools.combinations(range(16),10):
        a=H[:,cols];r=rank_mod(a);ranks[r]=ranks.get(r,0)+1
        if r<10:continue
        v=[((-1)**s)*int(round(np.linalg.det(np.delete(a,s,axis=0)))) for s in range(11)]
        gcd=math.gcd(*v);assert gcd
        v=[t//gcd for t in v]
        assert math.gcd(*v)==1
        assert all(sum(v[s]*int(a[s,j]) for s in range(11))==0 for j in range(10))
        dot=sum(t*(1<<s) for s,t in enumerate(v));assert dot
        val=(abs(dot)&-abs(dot)).bit_length()-1;hist[val]=hist.get(val,0)+1
        if val>maxv:maxv=val;example=dict(columns=cols,primitive_left_null=v,dot=dot)
    result=dict(status='EXACT_MODULAR_OBSTRUCTION',column_sets=math.comb(16,10),rank_histogram=ranks,
        two_adic_valuation_histogram=hist,max_valuation=maxv,example=example,
        conclusion='For denominators 2^n with n>=7, no <=10 exponent-Walsh terms represent 2^s/2^n modulo whole turns for s=0..10.',
        scope='separate Walsh selector per result bit; not a bound on general phase circuits or approximate selectors')
    Path(__file__).resolve().with_name('results').joinpath('modular_walsh_certificate.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
