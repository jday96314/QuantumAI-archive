"""Certify support search exactly over real-angle Walsh extensions.

The entries are +/-1. Eleven-column determinants have magnitude bounded by
11^(11/2)<2^20, and Cramer's numerators by 2048*11^5<2^29. A prime2^61-1
therefore detects singular matrices and zero coefficients exactly, without
floating-point rank decisions. Whole-turn lifts of the phase are not covered.
"""
import itertools,json,math
from pathlib import Path
P=(1<<61)-1
h=[[1-2*((s&m).bit_count()%2) for m in range(16)] for s in range(11)]

def solve(cols):
    a=[[h[s][m]%P for m in cols]+[(1<<s)%P] for s in range(11)]
    for j in range(11):
        piv=next((i for i in range(j,11) if a[i][j]),None)
        if piv is None:return None
        a[j],a[piv]=a[piv],a[j];iv=pow(a[j][j],-1,P)
        a[j][j:]=[t*iv%P for t in a[j][j:]]
        for i in range(j+1,11):
            f=a[i][j]
            if f:a[i][j:]=[(aa-f*bb)%P for aa,bb in zip(a[i][j:],a[j][j:])]
    v=[0]*11
    for i in reversed(range(11)):v[i]=(a[i][-1]-sum(a[i][j]*v[j] for j in range(i+1,11)))%P
    return v

if __name__=='__main__':
    hist={};sing=0
    for cols in itertools.combinations(range(16),11):
        v=solve(cols)
        if v is None:sing+=1
        else:
            n=sum(t!=0 for t in v);hist[n]=hist.get(n,0)+1
    result=dict(status='EXACT_CERTIFICATE',bases=math.comb(16,11),singular=sing,
        support_histogram=hist,minimum_support=min(hist),prime=P,
        scope='real-valued angle function on reachable s=0..10, unrestricted unreachable extension; excludes modular angle lifts')
    Path(__file__).resolve().with_name('results').joinpath('exact_walsh_support.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
