"""An unequal-width in-place adder using n-1 temporary ANDs.

Carries are computed from original operand bits, the sum is written, and
carries are erased by X-basis measurement while temporarily restoring each
original target bit. Missing high operand bits require no allocated zeros.
"""
import kickmix as km


def add(b,lhs,target,subtract=False,carry_in=None,carry_workspace=(),prefix_bits=0,qubit_limit=None):
    m,n=len(lhs),len(target)
    assert 0<m<=n
    if qubit_limit is not None:
        for candidate in range(min(m,n-1)+1):
            if candidate==0:extra=max(0,n-2-len(carry_workspace))
            else:
                initial=int(carry_in is None and not len(carry_workspace))
                remaining=len(carry_workspace)-int(carry_in is None and bool(len(carry_workspace)))
                extra=initial+max(0,n-candidate-2-remaining)
            if b.num_allocated_qubits+extra<=qubit_limit:
                prefix_bits=candidate;break
        else:raise ValueError(f'adder cannot fit {qubit_limit} qubits; already live {b.num_allocated_qubits}')
    k=min(prefix_bits,m,n-1)
    if k:
        if subtract:b.x(target)
        own=()
        available=carry_workspace
        if carry_in is not None:initial=carry_in
        elif len(available):initial=available[0];available=available[1:]
        else:own=b.alloc_qubits(1);initial=own[0]
        for i in range(k):
            c=initial if i==0 else lhs[i-1]
            b.cx(lhs[i],target[i]);b.cx(lhs[i],c);b.ccx(c,target[i],lhs[i])
        boundary=lhs[k-1]
        if m>k:add(b,lhs[k:],target[k:],carry_in=boundary,carry_workspace=available)
        else:add(b,km.array([boundary]),target[k:],carry_workspace=available)
        for i in reversed(range(k)):
            c=initial if i==0 else lhs[i-1]
            b.ccx(c,target[i],lhs[i]);b.cx(lhs[i],c);b.cx(c,target[i])
        if len(own):b.free(own)
        if subtract:b.x(target)
        return
    if subtract:b.x(target)
    borrowed=list(carry_workspace[:max(0,n-2)])
    extra_carries=b.alloc_qubits(max(0,n-2-len(borrowed)))
    carries=km.array([*borrowed,*extra_carries])
    for i in range(n-2):
        out=carries[i]
        if i==0 and carry_in is None:b.init_and(lhs[i],target[i],out)
        elif i>=m:b.init_and(carries[i-1],target[i],out)
        else:
            c=carry_in if i==0 else carries[i-1]
            b.cx(c,lhs[i]);b.cx(c,target[i])
            b.init_and(lhs[i],target[i],out)
            b.cx(c,out)
            b.cx(c,target[i]);b.cx(c,lhs[i])
    if n>=2:
        i=n-2
        if i==0 and carry_in is None:b.ccx(lhs[0],target[0],target[-1])
        elif i>=m:b.ccx(carries[i-1],target[i],target[-1])
        else:
            c=carry_in if i==0 else carries[i-1]
            b.cx(c,lhs[i]);b.cx(c,target[i])
            b.ccx(lhs[i],target[i],target[-1]);b.cx(c,target[-1])
            b.cx(c,target[i]);b.cx(c,lhs[i])
    for i in range(n):
        if i<m:b.cx(lhs[i],target[i])
        if i==0 and carry_in is not None:b.cx(carry_in,target[i])
        if 0<i<n-1:b.cx(carries[i-1],target[i])
    for i in reversed(range(n-2)):
        # Restore the original input bit for this carry's uncomputation.
        if i<m:b.cx(lhs[i],target[i])
        if i==0 and carry_in is not None:b.cx(carry_in,target[i])
        if i:b.cx(carries[i-1],target[i])
        out=carries[i]
        if i==0 and carry_in is None:b.del_and(lhs[i],target[i],out)
        elif i>=m:b.del_and(carries[i-1],target[i],out)
        else:
            c=carry_in if i==0 else carries[i-1]
            b.cx(c,out)
            b.cx(c,lhs[i]);b.cx(c,target[i])
            b.del_and(lhs[i],target[i],out)
            b.cx(c,target[i]);b.cx(c,lhs[i])
        if i<m:b.cx(lhs[i],target[i])
        if i==0 and carry_in is not None:b.cx(carry_in,target[i])
        if i:b.cx(carries[i-1],target[i])
    b.free(extra_carries)
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
