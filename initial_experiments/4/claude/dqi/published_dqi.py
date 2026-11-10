"""Runtime-mask adaptation of andigu/uncompute's Apache-2.0 DQI decoder.

See import_published.py and vendor/uncompute/LICENSE for upstream attribution.
No host verifier or reference data is changed.
"""
from pathlib import Path
import json, math, time
import numpy as np
import kickmix as km
from tower_inverse import map_ops,apply_map,linear_xor
from forney_dqi import product_terms
from import_published import load,HERE

def isomorphism(pub):
    native=km.GF2Field(12)
    def peval(x):
        v=0
        for i in reversed(range(13)):v=native.mul(v,x)^((pub.POLY>>i)&1)
        return v
    logs={native.pow(2,k):k for k in range(4095)}
    beta=next(x for x in range(2,4096) if peval(x)==0)
    k=logs[beta];r=pow(k,-1,4095)
    forward=[pub.gpow(pub.EXP[r],i) for i in range(12)]
    back=[native.pow(beta,i) for i in range(12)]
    def transform(v,cols):
        z=0
        for i in range(12):
            if v>>i&1:z^=cols[i]
        return z
    for a in range(4096):assert transform(transform(a,forward),back)==a
    return forward,back,transform

def sweep(b,pub,hook,masks,q,back,transform,positions):
    # Reuse public decoder's free wires when available; all are clean here.
    free=hook.get('free',[])
    own=[]
    def alloc(n):
        nonlocal free
        take=min(n,len(free));ids=free[-take:] if take else []
        if take:free=free[:-take]
        a=km.array([q[i] for i in ids])
        if take<n:
            extra=b.alloc_qubits(n-take);own.append(extra)
            a=km.array(list(a)+list(extra))
        return a
    root=alloc(12);omega=alloc(12);psi=alloc(12);u=alloc(1)[0]
    ladder=alloc(11)
    terms=product_terms(km.GF2Field(12))
    om=[[q[j] for j in w] for w in hook['omega']]
    eo=[[q[j] for j in w] for w in hook['eo']]
    ab=[[q[j] for j in w] for w in hook['ab']]
    def lin(w,target,c,square=False,native=False):
        cols=[pub.gmul(c,pub.gmul(1<<i,1<<i) if square else 1<<i) for i in range(12)]
        if native:cols=[transform(v,back) for v in cols]
        linear_xor(b,w,target,cols)
    for count,j in enumerate(positions):
        # Native position j corresponds to public alpha^(r*j); evaluate
        # at its reciprocal. The isomorphism's alpha image is forward[1].
        x=pub.gpow(pub.EXP[ISOM_R],(-j)%4095)
        def evals():
            for i,w in enumerate(om):lin(w,omega,pub.gpow(x,i),native=True)
            for i,w in enumerate(eo):
                h=i%18;c=pub.gpow(x,h)
                if i>=18:c=pub.gmul(c,pub.gsqrt(x))
                lin(w,root,c)
            for i,w in enumerate(ab):
                h=i%18;c=pub.gpow(x,2*h+(1 if i<18 else 0))
                lin(w,psi,c,square=True,native=True)
        evals()
        b.x(root);b.init_and(root[0],root[1],ladder[0])
        for k in range(1,11):b.init_and(ladder[k-1],root[k+1],ladder[k])
        def maskmap(jj):
            for k in range(12):
                sources=[omega[i] for i,l in terms[k] if l==jj]
                if not sources:continue
                b.push_condition(masks[j][k])
                b.cx(km.array(sources),u)
                b.pop_condition()
        for jj in range(12):
            maskmap(jj);b.ccz(ladder[10],u,psi[jj]);maskmap(jj)
        for k in reversed(range(1,11)):b.del_and(ladder[k-1],root[k+1],ladder[k])
        b.del_and(root[0],root[1],ladder[0]);b.x(root)
        evals()
        if count%512==0:print('sweep',count,flush=True)
    for reg in own:b.free(reg)

ISOM_R=None
def build_circuit(m=4095,n=70,chien=True,positions=None):
    global ISOM_R
    assert n==70
    pub=load();forward,back,transform=isomorphism(pub)
    ISOM_R=pub.LOG[forward[1]]
    hook=json.loads((HERE/'parameters/published_hook.json').read_text())
    rows=np.load(HERE/'parameters/published_tape.npy',mmap_mode='r')
    b=km.CircuitBuilder()
    syn=[b.create_quantum_register(12,f'S{i}') for i in range(70)]
    masks=[b.create_classical_register(12,f'E{i}') for i in range(m)]
    maxq=hook['maxq']
    if 'outcomes' not in hook:maxq+=1
    q=[v for w in syn for v in w]+list(b.alloc_qubits(maxq-840))
    fops=map_ops(forward)
    for w in syn:apply_map(b,w,fops)
    # Pool measurement bits at last use instead of importing monotonically
    # increasing outcome IDs. Outcome-only conditions remain free Clifford.
    sent=(1<<64)-1
    last={}
    for col,bit in [(4,8),(5,16)]:
        ids=np.flatnonzero((rows[:,0]>>(8+int(math.log2(bit))))&1==0)
        for i in ids:last[int(rows[i,col])]=max(int(i),last.get(int(rows[i,col]),-1))
    outs={};bits={}
    def outcome(i):
        if i not in outs:outs[i]=b.alloc_bits(1)[0]
        return outs[i]
    def cbit(i):
        if i not in bits:bits[i]=b.alloc_bits(1)[0]
        return bits[i]
    for index,row in enumerate(rows):
        if index==hook['index'] and chien:
            sweep(b,pub,hook,masks,q,back,transform,positions if positions is not None else range(m))
        tag=int(row[0]);kind=tag&255;mask=tag>>8
        v=[int(row[i+1]) if not(mask>>i&1) else sent for i in range(8)]
        a,c,t,ot,oc,rt,ct,cc=v
        if oc!=sent and kind!=10:b.push_condition(outcome(oc))
        if cc!=sent:b.push_condition(cbit(cc))
        if kind==0:b.x(q[t])
        elif kind==1:b.z(q[t])
        elif kind==2:b.neg()
        elif kind==3:b.cx(q[c],q[t])
        elif kind==4:b.cz(q[c],q[t])
        elif kind==5:b.swap(q[c],q[t])
        elif kind==6:b.ccx(q[a],q[c],q[t])
        elif kind==7:b.ccz(q[a],q[c],q[t])
        elif kind==8:b.reset(q[t])
        elif kind==9:b.hmr(q[t],outcome(ot))
        elif kind==10:b.cx(outcome(oc),outcome(ot))
        elif kind==11:b.x(cbit(ct))
        elif kind in [12,13]:b.bit_store(cbit(ct),kind==13)
        else:raise ValueError(kind)
        if cc!=sent:b.pop_condition()
        if oc!=sent and kind!=10:b.pop_condition()
        for i in set([ot,oc])-{sent}:
            if last[i]==index:
                b.bit_store(outs[i],False);b.free(outs.pop(i))
        if index%1000000==0:print('translate',index,flush=True)
    for w in syn:apply_map(b,w,fops,True)
    return b.finish_circuit()

if __name__=='__main__':
    import sys
    start=time.time();c=build_circuit()
    print('built',c.num_qubits,c.operation_statistics().num_toffoli_like,time.time()-start,flush=True)
    c.to_file(sys.argv[1],format='kmb')
