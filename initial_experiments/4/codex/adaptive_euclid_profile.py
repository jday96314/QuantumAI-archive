"""Exact stage counts and width-constrained schedules for the compact core."""
import json
from pathlib import Path
from import_published import load
from blocked_euclid_sequence import frames
from blocked_monic_euclid import apply as block
from monomial_serial_euclid import apply as serial
from monic_euclid_local import first_step
from deferred_euclid_padding import direct_repad
from mobius_series import BoundedArithmetic


def profile(N,k=None,width=None,finish=False):
    m=load();b=m.Builder();words=[b.qubits(12) for _ in range(N)];flags=b.qubits(N//2)
    class Counter(m.Em):
        T=0
        def emit_mac(self,rt):self.T+=sum(t[0]=='P' for t in rt)
        def emit_free(self,rt):pass
    E=Counter(b,*b.qubits(3),m.Alloc(b));math=BoundedArithmetic(m,E)
    if k is None:
        if finish:direct_repad(m,E,words,flags)
        else:
            with m.eq_test(E,words[-1],0) as zero:E.cx(zero,flags[0])
            E.cx(flags[0],words[-1][0]);first_step(m,E,math,words)
    else:
        stages,_=frames(words,1);_,_,bp,bb,c,p=stages[k-1]
        if width==1:serial(m,E,math,bp,bb,c,p,flags[k])
        else:block(m,E,math,bp,bb,c,p,width,flags=flags[k:k+width],monomial=True,suffix=True)
    return dict(k=k,width=width,Q=b.next_qubit,T=E.T+b._cur_stats['toffoli'])


def main():
    N=70;rows=[]
    for k in range(1,N//2):
        for width in [1,2,3,4,5,6]:
            if k+width>N//2:continue
            r=profile(N,k,width);rows.append(r)
        print(json.dumps(dict(k=k,rows=rows[-min(6,N//2-k):])),flush=True)
    first=profile(N);last=profile(N,finish=True);schedules=[]
    for cap in sorted({r['Q'] for r in rows}):
        dp={N//2:(0,[])}
        for k in reversed(range(1,N//2)):
            choices=[(r['T']+dp[k+r['width']][0],[r]+dp[k+r['width']][1]) for r in rows
                     if r['k']==k and r['Q']<=cap and k+r['width'] in dp]
            if choices:dp[k]=min(choices,key=lambda x:x[0])
        if 1 in dp:
            cost,path=dp[1];r=dict(Q=max(cap,first['Q'],last['Q']),T=cost+first['T']+last['T'],stages=path);schedules.append(r)
    out=dict(scope='compact forward core including monomial completion and direct repadding; excludes scan, repetition and full decoder validation',N=N,first=first,last=last,rows=rows,schedules=schedules)
    Path('results/adaptive_euclid_profiles_n70.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(schedules=schedules)),flush=True)
if __name__=='__main__':main()
