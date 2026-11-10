"""Recover the phase from the retained last 35 syndromes and the locator."""
import kickmix as km
from fast_dqi import FastBuilder
from noanc_multiplier import NoAncMultiply
from recursive_forney import RecursiveForney
from tower_inverse import apply_map,map_ops
from toom_dqi import plan,lin_xor

class Kernel(NoAncMultiply,RecursiveForney):
    gf=km.GF2Field(12)
    t=35

K=Kernel();GF=K.gf

def middle(b,sigma,low,high):
    """Transpose a Toom product: high[i] ^= sum_j sigma[L+i-j] low[j]."""
    n=len(low)+len(high)-1
    _,E,I=plan(n)
    s=b.alloc_qubits(12);x=b.alloc_qubits(12);v=b.alloc_qubits(12)
    for a in range(n):
        se=[(I[a][j],sigma[j+1],s) for j in range(n)]
        xe=[(E[a][j],low[-1-j],x) for j in range(len(low))]
        lin_xor(b,se);lin_xor(b,xe);K._mul(b,s,x,v)
        lin_xor(b,[(E[a][i],v,high[i]) for i in range(len(high))])
        b.del_gf2_mul(s,x,target=v,field=GF)
        lin_xor(b,xe);lin_xor(b,se)
    for reg in [s,x,v]:b.free(reg)

def omega(b,sigma,syn,inverse=False):
    if len(syn)<=1:return
    mid=len(syn)//2;lo=syn[:mid];hi=syn[mid:]
    if inverse:
        omega(b,sigma,lo,True);middle(b,sigma,lo,hi);omega(b,sigma,hi,True)
    else:
        omega(b,sigma,hi);middle(b,sigma,lo,hi);omega(b,sigma,lo)

def prefix_zero(b,words,target):
    if len(words)>8:
        chunks=[words[i:i+8] for i in range(0,len(words),8)]
        flags=b.alloc_qubits(len(chunks))
        for chunk,flag in zip(chunks,flags):prefix_zero(b,chunk,flag)
        b.flip_if_equal(flags,(1<<len(chunks))-1,target=target)
        for chunk,flag in reversed(list(zip(chunks,flags))):prefix_zero(b,chunk,flag)
        b.free(flags)
        return
    flags=b.alloc_qubits(len(words))
    for w,f in zip(words,flags):b.flip_if_equal(w,0,target=f)
    b.flip_if_equal(flags,(1<<len(words))-1,target=target)
    for w,f in reversed(list(zip(words,flags))):b.flip_if_equal(w,0,target=f)
    b.free(flags)

def word_mul(b,x,scale,eraser):
    v=b.alloc_qubits(12);K._mul(b,x,scale,v)
    mask=b.alloc_bits(12);b.hmr(x,mask)
    b.gf2_phase_by_product(v,eraser,mask=mask,field=GF)
    b.bit_store(mask,False);b.free(mask);b.swap(x,v);b.free(v)

def chien(b,sigma,om,masks):
    # The reversed syndrome sequence has reciprocal error locations.
    step=GF.primitive_element
    # The normalized constant coefficient is implicitly one. Keep its original
    # scale in sigma[0] and use sigma[2] to sum the nonconstant terms instead.
    acc_d,acc_s,acc_o=sigma[1],sigma[2],om[0]
    odd=[sigma[d] for d in range(3,36,2)];even=[sigma[d] for d in range(4,36,2)]
    src=km.array([q for w in odd+[sigma[1]]+even+om[1:] for q in w])
    dst=km.array([q for w in [acc_d]*len(odd)+[acc_s]*(1+len(even))+[acc_o]*34 for q in w])
    rsrc=km.array(list(reversed(list(src))));rdst=km.array(list(reversed(list(dst))))
    for j,mask in enumerate(masks):
        point=GF.pow(step,j)
        b.cx(src,dst)
        if point!=1:b.gf2_imul(point,target=acc_o,field=GF)
        b.x(acc_s[0])
        K._forney_kernel(b,acc_s,acc_d,acc_o,mask)
        b.x(acc_s[0])
        if point!=1:b.gf2_idiv(point,target=acc_o,field=GF)
        b.cx(rsrc,rdst)
        if j+1<len(masks):
            for i in range(1,36):b.gf2_imul(GF.pow(step,i),target=sigma[i],field=GF)
            for i in range(1,35):b.gf2_imul(GF.pow(step,i),target=om[i],field=GF)
    if len(masks)>1:
        for i in range(1,36):b.gf2_idiv(GF.pow(step,i*(len(masks)-1)%4095),target=sigma[i],field=GF)
        for i in range(1,35):b.gf2_idiv(GF.pow(step,i*(len(masks)-1)%4095),target=om[i],field=GF)

def phase_kernel(builder,hook,masks,q,forward,back):
    b=FastBuilder(builder)
    sigma=[km.array([q[j] for j in w]) for w in hook['sigma']]
    syn=[km.array([q[j] for j in w]) for w in hook['tail']]
    backops=map_ops(back)
    for w in sigma+syn:apply_map(b,w,backops)
    sigma=list(reversed(sigma));flags=b.alloc_qubits(6)
    for k in reversed(range(6)):
        shift=1<<k
        prefix_zero(b,sigma[:shift],flags[k])
        for i in range(36-shift):b.cswap(flags[k],sigma[i],sigma[i+shift])
    inv,sc=b.init_gf2_inverse_with_scaffold(sigma[0],field=GF)
    for w in sigma[1:]:word_mul(b,w,inv,sigma[0])
    b.del_gf2_inverse_with_scaffold(sigma[0],target=inv,scaffold=sc,field=GF);b.free(inv);b.free(sc)
    omega(b,sigma,syn)
    maskops=[]
    for j,mask in enumerate(masks):
        c=GF.pow(GF.primitive_element,(-71*j)%4095)
        columns=tuple(sum(((GF.mul(c,1<<i)>>k)&1)<<i for i in range(12)) for k in range(12))
        ops=map_ops(columns);maskops.append(ops);apply_map(b,mask,ops)
    chien(b,sigma,syn,masks)
    for mask,ops in zip(masks,maskops):apply_map(b,mask,ops,True)
    omega(b,sigma,syn,True)
    inv,sc=b.init_gf2_inverse_with_scaffold(sigma[0],field=GF)
    for w in reversed(sigma[1:]):word_mul(b,w,sigma[0],inv)
    b.del_gf2_inverse_with_scaffold(sigma[0],target=inv,scaffold=sc,field=GF);b.free(inv);b.free(sc)
    for k in range(6):
        shift=1<<k
        for i in reversed(range(36-shift)):b.cswap(flags[k],sigma[i],sigma[i+shift])
        prefix_zero(b,sigma[:shift],flags[k])
    b.free(flags)
    for w in list(reversed(sigma))+syn:apply_map(b,w,backops,True)

def test_omega():
    import numpy as np
    rng=np.random.default_rng(50)
    for n in [2,3,6,17,35]:
        b=FastBuilder(km.CircuitBuilder());sig=[b.create_quantum_register(12,f's{i}') for i in range(n)]
        x=[b.create_quantum_register(12,f'x{i}') for i in range(n)]
        omega(b,sig,x);c=b.finish_circuit();N=256
        ss=rng.integers(0,4096,(n,N));ss[0]=1;xx=rng.integers(0,4096,(n,N))
        sim=km.Simulator(batch_size=N,seed=30);sim.use_same_registers_as(c)
        bits=lambda a:((a[None,:]>>np.arange(12)[:,None])&1).astype(bool)
        for w,v in zip(sig,ss):sim.write_across_shots(w,bits(v))
        for w,v in zip(x,xx):sim.write_across_shots(w,bits(v))
        sim.do(c)
        for i in range(n):
            expected=[]
            for s in range(N):
                z=0
                for j in range(i+1):z^=GF.mul(int(ss[j,s]),int(xx[i-j,s]))
                expected.append(z)
            assert np.array_equal(sim.read_across_shots(x[i]),bits(np.array(expected)))
        assert not sim.read_phase_flipped_across_shots().any()
        print('omega PASS',n,c.num_qubits,c.operation_statistics().num_toffoli_like,flush=True)

if __name__=='__main__':test_omega()
