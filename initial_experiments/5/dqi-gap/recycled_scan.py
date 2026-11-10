"""Reversibly clear unary completion flags for use as scan accumulators.

For a non-unary residual path these accumulators can be dirty. Evaluation and
the phase kernel restore their inputs exactly, then the XOR decoder restores
the original flags. Thus approximate decoding never destroys history.
"""
import kickmix as km


def toggle_unary(b,flags,counter):
    n=len(flags);ladder=b.alloc_qubits(len(counter)-1)
    for value in range(1,n+1):
        for i,q in enumerate(counter):
            if not(value>>i&1):b.x(q)
        b.init_and(counter[0],counter[1],ladder[0])
        for i in range(2,len(counter)):b.init_and(ladder[i-2],counter[i],ladder[i-1])
        for flag in flags[n-value:]:b.cx(ladder[-1],flag)
        for i in reversed(range(2,len(counter))):b.del_and(ladder[i-2],counter[i],ladder[i-1])
        b.del_and(counter[0],counter[1],ladder[0])
        for i,q in enumerate(counter):
            if not(value>>i&1):b.x(q)
    b.free(ladder)


def test():
    import json,numpy as np
    from fast_dqi import FastBuilder
    b=FastBuilder(km.CircuitBuilder());flags=b.create_quantum_register(35,'flags');counter=b.create_quantum_register(6,'counter')
    toggle_unary(b,flags,counter);c=b.finish_circuit();s=km.Simulator(batch_size=256,seed=51551);s.use_same_registers_as(c)
    rng=np.random.default_rng(91524);cv=np.arange(256)%64;fv=rng.integers(0,2,(35,256),dtype=bool)
    fv[:,::2]=np.arange(35)[:,None]>=35-cv[None,::2]
    s.write_across_shots(flags,fv);s.write_across_shots(counter,((cv[None,:]>>np.arange(6)[:,None])&1).astype(bool));s.do(c)
    want=fv.copy()
    for i,v in enumerate(cv):
        if v<=35:want[35-v:,i]^=True
    assert np.array_equal(s.read_across_shots(flags),want)
    assert not s.read_phase_flipped_across_shots().any();s.do(c)
    assert np.array_equal(s.read_across_shots(flags),fv);assert not s.read_phase_flipped_across_shots().any()
    print(json.dumps(dict(status='PASS',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,shots=256)),flush=True)
if __name__=='__main__':test()
