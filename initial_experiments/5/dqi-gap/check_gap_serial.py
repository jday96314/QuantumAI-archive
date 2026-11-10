"""Independent polynomial and native-field tests of all serial update charts."""
import argparse,json
from pathlib import Path
import numpy as np
from import_published import load
from low_workspace_inverse import bind
from mobius_series import BoundedArithmetic
from gap_serial import apply
from native_tools import finish,check
from classical_repair import mul,inv,product,plus
from check_gap_block import sp,sd
from published_dqi import isomorphism


def reference(bp,b,c,p,stash,carry,wait,shift):
    n=len(c);flags=0;old=list(b);oldstash=list(stash)
    if wait:
        ns=list(reversed(c))[:2];np_=[0]+bp;nb=[0]+b
    elif carry:
        f=list(stash)+list(reversed(c))[:2]
        f=sd(sp(f,b[::-1],4),[inv(p),0]+sp(f,bp[::-1],2),4)
        d=f[1];q1=mul(f[2],inv(d));q0=mul(f[3],inv(d))^mul(q1,q1)
        nb=plus(product([q0,q1,1],b[1:]),[mul(d,x) for x in bp[1:]] if len(bp)>1 else [])
        np_=list(b);p=mul(p,inv(d));ns=[0,0]
    else:
        f=list(reversed(c))[:2];f=sd(sp(f,b[::-1],2),[inv(p),0],2)
        flags=int(f[0]==0);beta=f[0] or 1;alpha=mul(f[1],inv(beta))
        if flags:
            k=len(b)-1
            fake=plus(product([alpha,1],b),bp)
            fake += [0]*max(0,k+2-len(fake))
            for i in range(1,k-1):fake[i]^=b[i-1]
            if k>1:fake[k-1]^=b[k-2]
            fake[k]^=b[k-1];fake[k-1]^=1
            np_=[fake[k]]+fake[:k-1]+[1]
            nb=[fake[k-1]]+b[:-1]+[1]
        else:
            np_=list(b);nb=plus(product([alpha,1],b),[mul(beta,x) for x in (bp[1:] if shift else bp)])
        p=mul(p,inv(beta));ns=list(stash)
    numerator=sp(nb[::-1],c[::-1],n)[2:]
    if carry:
        corr=sp(nb[::-1],oldstash,n+2)[4:]
        numerator=[a^b for a,b in zip(numerator,corr)]
    if wait:
        corr=sp(old[::-1],ns,n)[2:]
        numerator=[a^b for a,b in zip(numerator,corr)]
    nc=sd(numerator,old[::-1],n-2)[::-1]
    return np_,nb,nc,p,ns,flags


def build(N,k,which):
    m=load();m.inv_word=bind(m);b=m.Builder();W=[b.qubits(12) for _ in range(N)];Z=[b.qubits(12) for _ in range(2)];F=b.qubits(4)
    E=m.Em(b,*b.qubits(3),m.Alloc(b));math=BoundedArithmetic(m,E);math.virtual_zero=True
    bp=W[:k-1];bb=W[k-1:2*k-1];c=W[2*k-1:-1];p=W[-1]
    for back in ([False,True] if which=='both' else [True,False] if which=='inverse_both' else [which=='inverse']):
        apply(m,E,math,bp,bb,c,p,Z,*F,back)
    regs=[(f'W{i}',w) for i,w in enumerate(W)]+[(f'Z{i}',w) for i,w in enumerate(Z)]+[('F',F)]
    return finish(b,regs,f'gap_serial_n{N}_k{k}_{which}')


def test(N,k,shots):
    rng=np.random.default_rng(913531+N+k);m=load();fw,bw,tr=isomorphism(m);starts=[];wants=[];groups=[]
    while len(starts)<shots:
        row=rng.integers(0,4096,N).tolist();row[-1]=int(rng.integers(1,4096));stash=[0,0]
        kind=len(starts)%(6 if k>1 else 4)
        carry=int(kind==1);wait=int(kind==2);shift=int(kind in [4,5]);completion=kind in [3,5]
        if carry:
            if k>1:row[0]=0
            row[k-1]=0;stash[1]=int(rng.integers(1,4096))
            if k==1:row[-1]=1
        elif wait:row[-2]=0;row[-3]=int(rng.integers(1,4096))
        elif completion:row[2*k-1:-1]=[0]*(N-2*k)
        else:row[-2]=int(rng.integers(1,4096))
        if shift:row[0]=0
        bp=row[:k-1]+[1];bb=row[k-1:2*k-1]+[1];cc=row[2*k-1:-1];p=row[-1]
        np_,nb,nc,npi,ns,flag=reference(bp,bb,cc,p,stash,carry,wait,shift)
        out=row.copy();newpi=list(range(k-1,2*k-1));newbi=list(range(k-1))+[N-2,N-3]
        assert len(np_)==k+1 and len(nb)==k+2,(k,kind,np_,nb)
        for i,v in zip(newpi,np_[:-1]):out[i]=v
        for i,v in zip(newbi,nb[:-1]):out[i]=v
        out[2*k-1:-3]=nc;out[-1]=npi
        starts.append(row+stash+[carry+2*wait+4*shift]);wants.append(out+ns+[carry+2*wait+4*shift+8*flag]);groups.append(kind)
    convert=lambda rows:{**{f'W{i}':np.array([tr(int(r[i]),fw) for r in rows]) for i in range(N)},**{f'Z{i}':np.array([tr(int(r[N+i]),fw) for r in rows]) for i in range(2)},'F':np.array([r[-1] for r in rows])}
    ss,ww=convert(starts),convert(wants);reports=[]
    for which in ['forward','inverse','both','inverse_both']:
        c=build(N,k,which);src=ww if which in ['inverse','inverse_both'] else ss;dst=ww if which in ['forward','inverse_both'] else ss
        report=check(c,src,dst,shots);report['which']=which;reports.append(report);print(json.dumps(report),flush=True)
    result=dict(N=N,k=k,shots=shots,status='PASS' if all(r['status']=='PASS' for r in reports) else 'FAIL',groups=groups,circuits=reports)
    Path(f'results/gap_serial_n{N}_k{k}.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=14);p.add_argument('--k',type=int,default=3);p.add_argument('--shots',type=int,default=256);p.add_argument('--unified',action='store_true');p.add_argument('--onepass',action='store_true');a=p.parse_args()
    if a.unified:
        import gap_serial
        gap_serial.UNIFIED=True
    if a.onepass:
        from onepass_gap import apply
    test(a.n,a.k,a.shots)
