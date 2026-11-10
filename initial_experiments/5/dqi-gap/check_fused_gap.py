"""Independent full-state tests for the unified block, including boundary carry."""
import argparse,json
from pathlib import Path
import numpy as np
from import_published import load
from low_workspace_inverse import bind
from mobius_series import BoundedArithmetic
from blocked_monic_euclid import schedule
from native_tools import finish,check
from classical_repair import mul,inv
from check_gap_block import sp,sd
from check_gap_serial import reference
from fused_gap_block import apply
from published_dqi import isomorphism


def build(N,k,count,which):
    m=load();m.inv_word=bind(m);b=m.Builder();W=[b.qubits(12) for _ in range(N)];Z=[b.qubits(12) for _ in range(2)];F=b.qubits(2*(k+count)+1)
    waits=F[:k+count];flags=F[k+count:-1];terminal=F[-1]
    E=m.Em(b,*b.qubits(3),m.Alloc(b));math=BoundedArithmetic(m,E);math.virtual_zero=True
    bp=W[:k-1];bb=W[k-1:2*k-1];c=W[2*k-1:-1];p=W[-1]
    for back in ([False,True] if which=='both' else [True,False] if which=='inverse_both' else [which=='inverse']):
        apply(m,E,math,bp,bb,c,p,Z,k,count,waits,flags,terminal,back)
    regs=[(f'W{i}',w) for i,w in enumerate(W)]+[(f'Z{i}',w) for i,w in enumerate(Z)]+[('F',F)]
    return finish(b,regs,f'fused_gap_n{N}_k{k}_b{count}_{which}')


def test(N,k,count,shots):
    rng=np.random.default_rng(884051+N+k);m=load();fw,bw,tr=isomorphism(m);starts=[];wants=[];groups=[]
    bpidx=list(range(k-1));bidx=list(range(k-1,2*k-1));cidx=list(range(2*k-1,N-1));rawidx=cidx[::-1][:2*count]
    _,(newpi,newbi)=schedule(bpidx,bidx,rawidx)
    while len(starts)<shots:
        kinds=[0]+list(range(1,count+1))+[4,5,6,7]
        if k==1:kinds=[v for v in kinds if v not in [5,7]]
        row=rng.integers(0,4096,N).tolist();p=row[-1]=int(rng.integers(1,4096));kind=kinds[len(starts)%len(kinds)]
        carry=int(kind==4);shift=int(kind in [5,7]);completion=kind in [6,7]
        if carry:
            row[0]=row[k-1]=0
            if k==1:p=row[-1]=1
        if shift:row[0]=0
        bp=row[:k-1]+[1];bb=row[k-1:2*k-1]+[1]
        nn=2*count+(2 if carry else 0);f=rng.integers(0,4096,nn).tolist();f[0]=int(rng.integers(1,4096))
        if kind in [1,2,3]:
            gap=kind-1;f=rng.integers(0,4096,2*(count-gap)).tolist();f[0]=0;f[1]=int(rng.integers(1,4096))
            for _ in range(gap):
                be=int(rng.integers(1,4096));al=int(rng.integers(0,4096));f=sd([be],[1,al]+f,len(f)+2)
        if carry:f[0]=0;f[1]=int(rng.integers(1,4096))
        skip=3 if shift else 2;den=bb[::-1]+[0]*max(0,nn-len(bb))
        add=sp(f,[] if carry and k==1 else bp[::-1],max(0,nn-skip))
        for i,v in enumerate(add):den[i+skip]^=v
        t=sd([mul(inv(p),v) for v in f],den,nn)
        stash=t[:2] if carry else [0,0];head=t[2:] if carry else t
        cc=row[2*k-1:-1];cc[-2*count:]=head[::-1]
        if completion:cc=[0]*len(cc)
        row[2*k-1:-1]=cc
        waits=[0]*(k+count);flags=[0]*(k+count);waits[k-1]=carry
        if k>=2:waits[k-2]=shift
        term=0
        initial_flags=sum(v<<i for i,v in enumerate(waits+flags+[term]))
        np_,nb,nc,npi,ns=bp,bb,cc,p,stash
        for j in range(count):
            kk=k+j;ca=waits[kk-1];sh=(waits[kk-2] if kk>=2 else 0)^term
            wa=int(not(nc[-1] or not nc[-2] or ca or sh))
            np_,nb,nc,npi,ns,flag=reference(np_,nb,nc,npi,ns,ca,wa,sh)
            waits[kk]=wa;flags[kk]=flag;term^=(waits[kk-2] if kk>=2 else 0)&flag
        # The block returns the local completion/wait flags. Its caller
        # updates the persistent terminal-shift bit once per block.
        final_flags=sum(v<<i for i,v in enumerate(waits+flags+[0]))
        out=row.copy()
        for i,v in zip(newpi,np_[:-1]):out[i]=v
        for i,v in zip(newbi,nb[:-1]):out[i]=v
        for i,v in zip(cidx[:-2*count],nc):out[i]=v
        out[-1]=npi;starts.append(row+stash+[initial_flags]);wants.append(out+ns+[final_flags]);groups.append(kind)
    convert=lambda rows:{**{f'W{i}':np.array([tr(int(r[i]),fw) for r in rows]) for i in range(N)},**{f'Z{i}':np.array([tr(int(r[N+i]),fw) for r in rows]) for i in range(2)},'F':np.array([r[-1] for r in rows])}
    ss,ww=convert(starts),convert(wants);reports=[]
    for which in ['forward','inverse','both','inverse_both']:
        c=build(N,k,count,which);src=ww if which in ['inverse','inverse_both'] else ss;dst=ww if which in ['forward','inverse_both'] else ss
        r=check(c,src,dst,shots);r['which']=which;reports.append(r);print(json.dumps(r),flush=True)
    report=dict(N=N,k=k,count=count,shots=shots,status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',groups=groups,circuits=reports)
    Path(f'results/fused_gap_n{N}_k{k}_b{count}.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=22);p.add_argument('--k',type=int,default=3);p.add_argument('--count',type=int,default=3);p.add_argument('--shots',type=int,default=128);a=p.parse_args();test(a.n,a.k,a.count,a.shots)
