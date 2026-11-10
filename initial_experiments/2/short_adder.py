"""An unequal-width in-place adder using n-1 temporary ANDs.

Carries are computed from original operand bits, the sum is written, and
carries are erased by X-basis measurement while temporarily restoring each
original target bit. Missing high operand bits require no allocated zeros.
"""
import kickmix as km


def add(b,lhs,target,subtract=False):
    m,n=len(lhs),len(target)
    assert 0<m<=n
    if subtract:b.x(target)
    carries=b.alloc_qubits(max(0,n-2))
    for i in range(n-2):
        out=carries[i]
        if i==0:b.init_and(lhs[i],target[i],out)
        elif i>=m:b.init_and(carries[i-1],target[i],out)
        else:
            c=carries[i-1]
            b.cx(c,lhs[i]);b.cx(c,target[i])
            b.init_and(lhs[i],target[i],out)
            b.cx(c,out)
            b.cx(c,target[i]);b.cx(c,lhs[i])
    if n>=2:
        i=n-2
        if i==0:b.ccx(lhs[0],target[0],target[-1])
        elif i>=m:b.ccx(carries[i-1],target[i],target[-1])
        else:
            c=carries[i-1]
            b.cx(c,lhs[i]);b.cx(c,target[i])
            b.ccx(lhs[i],target[i],target[-1]);b.cx(c,target[-1])
            b.cx(c,target[i]);b.cx(c,lhs[i])
    for i in range(n):
        if i<m:b.cx(lhs[i],target[i])
        if 0<i<n-1:b.cx(carries[i-1],target[i])
    for i in reversed(range(n-2)):
        # Restore the original input bit for this carry's uncomputation.
        if i<m:b.cx(lhs[i],target[i])
        if i:b.cx(carries[i-1],target[i])
        out=carries[i]
        if i==0:b.del_and(lhs[i],target[i],out)
        elif i>=m:b.del_and(carries[i-1],target[i],out)
        else:
            c=carries[i-1]
            b.cx(c,out)
            b.cx(c,lhs[i]);b.cx(c,target[i])
            b.del_and(lhs[i],target[i],out)
            b.cx(c,target[i]);b.cx(c,lhs[i])
        if i<m:b.cx(lhs[i],target[i])
        if i:b.cx(carries[i-1],target[i])
    b.free(carries)
    if subtract:b.x(target)


def test():
    import numpy as np
    rng=np.random.default_rng(20261022)
    for m,n in [(1,1),(1,6),(3,6),(6,6),(20,28)]:
        for sub in [False,True]:
            b=km.CircuitBuilder();a=b.create_quantum_register(m,'a');t=b.create_quantum_register(n,'t')
            add(b,a,t,sub);c=b.finish_circuit();N=4096
            x=rng.integers(0,1<<m,N);y=rng.integers(0,1<<n,N)
            if n<=6:
                pairs=np.array([(i,j) for i in range(1<<m) for j in range(1<<n)])
                x=np.resize(pairs[:,0],N);y=np.resize(pairs[:,1],N)
            bits=lambda z,w:((z[None,:]>>np.arange(w)[:,None])&1).astype(bool)
            s=km.Simulator(batch_size=N,seed=17);s.use_same_registers_as(c)
            s.write_across_shots(a,bits(x,m));s.write_across_shots(t,bits(y,n));s.do(c)
            expected=(y-x if sub else y+x) % (1<<n)
            assert np.array_equal(s.read_across_shots(a),bits(x,m)),(m,n,sub)
            assert np.array_equal(s.read_across_shots(t),bits(expected,n)),(m,n,sub)
            assert not s.read_phase_flipped_across_shots().any(),(m,n,sub)
            assert not s.read_across_shots(km.array([km.q(i) for i in range(m+n,c.num_qubits)])).any(),(m,n,sub)
            print('PASS',m,n,'subtract',sub,'Toffolis',c.operation_statistics().num_toffoli_like,flush=True)
if __name__=='__main__':test()
