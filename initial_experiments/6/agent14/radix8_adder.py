"""Radix-eight signed multiplication into a modular accumulator.

Two adjacent multiplier bits share one ripple addition. Partial products use
measurement erasure; the multiplier and multiplicand are restored exactly.
"""
import kickmix as km
from short_adder import add

def accumulate(b,source,coordinate,target,exponent,*,subtract=False,rounded=True,borrow=(),auxiliary=(),adder_prefix_bits=0,qubit_limit=None,product_limit=None):
    L=len(coordinate);pad=b.alloc_qubits((-L)%3)
    for q in pad:b.cx(coordinate[-1],q)
    coord=km.array([*coordinate,*pad])
    for g in range(len(coord)//3):
        e=exponent+3*g;shapes=[]
        for branch in range(3):
            ee=e+branch
            shape=[(i,i+ee) for i in range(len(source)) if 0<=i+ee<len(target)]
            rs=-ee-1 if rounded and ee<0 and -ee-1<len(source) else None
            shapes.append((shape,rs))
        if not any(shape or rs is not None for shape,rs in shapes):continue
        current,middle,sign=coord[3*g:3*g+3];prev=coord[3*g-1] if g else None
        extra_tmp=b.alloc_qubits(max(0,2-len(auxiliary)));tmp=km.array([*auxiliary[:2],*extra_tmp]);diff,mag4=tmp
        def make_a():
            # a=(prev xor sign)*(original current xor sign).
            b.cx(sign,current)
            if prev is None:b.init_and(sign,current,diff)
            else:
                b.cx(sign,prev);b.init_and(prev,current,diff);b.cx(sign,prev)
            b.cx(sign,current)
        def erase_a():
            b.cx(sign,current)
            if prev is None:b.del_and(sign,current,diff)
            else:
                b.cx(sign,prev);b.del_and(prev,current,diff);b.cx(sign,prev)
            b.cx(sign,current)
        make_a();b.cx(sign,middle);b.init_and(diff,middle,mag4);b.cx(sign,middle);erase_a()
        if prev is not None:b.cx(prev,current);b.cx(prev,middle)
        b.cx(sign,diff)
        if prev is not None:b.cx(prev,diff)
        b.ccx(current,diff,middle)
        if prev is not None:b.cx(prev,diff)
        b.cx(sign,diff)
        controls=[current,middle,mag4];top=False
        normal={};carry=[]
        for ctrl,(shape,rs) in zip(controls,shapes):
            for i,pos in shape:normal.setdefault(pos,[]).append((ctrl,source[i]))
            if rs is not None:carry.append((ctrl,source[rs]))
        if normal:
            low,high=min(normal),max(normal);count=high-low+1
        else:low=0;count=0
        full_count=count
        if product_limit is not None:count=min(count,product_limit)
        total=count+bool(carry)
        borrowed=list(borrow[:min(total,len(borrow))]);extra=b.alloc_qubits(total-len(borrowed))
        products=km.array([*borrowed,*extra]);claims=[]
        for k in range(count):claims.append(normal.get(low+k,[]))
        if carry:claims.append(carry)
        for out,terms in zip(products,claims):
            for ctrl,src in terms:b.ccx(ctrl,src,out)
        lane=target[low:]
        b.cx(sign,lane)
        carry_work=km.array([*([diff] if not top else []),*borrow[len(borrowed):]])
        if count:add(b,products[:count],lane,subtract=subtract,carry_in=products[-1] if carry else None,carry_workspace=carry_work,prefix_bits=adder_prefix_bits,qubit_limit=qubit_limit)
        else:add(b,products,lane,subtract=subtract,carry_workspace=carry_work,prefix_bits=adder_prefix_bits,qubit_limit=qubit_limit)
        meas=b.alloc_bits(total);b.hmr(products,meas)
        for bit,terms in zip(meas,claims):
            b.push_condition(bit)
            for ctrl,src in terms:b.cz(ctrl,src)
            b.pop_condition()
        b.bit_store(meas,False);b.free(meas);b.free(extra)
        # The high product bits can be added after erasing the low products.
        # The two Booth magnitude controls are disjoint, so their XOR equals
        # their integer sum. Each deferred term needs only one clean bit.
        for k in range(count,full_count):
            for ctrl,src in normal.get(low+k,[]):
                p=b.alloc_qubits(1)
                b.init_and(ctrl,src,p[0])
                add(b,p,target[low+k:],subtract=subtract,qubit_limit=qubit_limit)
                b.del_and(ctrl,src,p[0]);b.free(p)
        b.cx(sign,lane)
        b.cx(sign,diff)
        if prev is not None:b.cx(prev,diff)
        b.ccx(current,diff,middle)
        if prev is not None:b.cx(prev,diff)
        b.cx(sign,diff)
        if prev is not None:b.cx(prev,middle);b.cx(prev,current)
        make_a();b.cx(sign,middle);b.del_and(diff,middle,mag4);b.cx(sign,middle);erase_a()
        b.free(extra_tmp)
    for q in pad:b.cx(coordinate[-1],q)
    b.free(pad)

def test():
    import numpy as np
    rng=np.random.default_rng(8102)
    for ns,L,nt,ex in [(5,5,8,-2),(5,6,9,-4),(4,5,7,0),(3,4,6,-6)]:
        for sub in [False,True]:
            b=km.CircuitBuilder();src=b.create_quantum_register(ns,'s');z=b.create_quantum_register(L,'z');t=b.create_quantum_register(nt,'t')
            accumulate(b,src,z,t,ex,subtract=sub);c=b.finish_circuit();N=4096
            aa=np.resize(np.repeat(np.arange(1<<ns),1<<L),N);zz=np.resize(np.tile(np.arange(1<<L),1<<ns),N);tt=rng.integers(0,1<<nt,N)
            delta=np.zeros(N,dtype=np.int64)
            for g in range((L+2)//3):
                bit=lambda k:0 if k<0 else (zz>>min(k,L-1))&1
                digit=bit(3*g-1)+bit(3*g)+2*bit(3*g+1)-4*bit(3*g+2)
                for mag in [1,2,4]:
                    e=ex+3*g+(mag.bit_length()-1);v=aa<<e if e>=0 else (aa>>-e)+((aa>>(-e-1))&1)
                    delta+=np.where((abs(digit)&mag)!=0,np.sign(digit)*v,0)
            expected=(tt-delta if sub else tt+delta)%(1<<nt)
            bits=lambda x,w:((x[None,:]>>np.arange(w)[:,None])&1).astype(bool)
            sim=km.Simulator(batch_size=N,seed=181);sim.use_same_registers_as(c)
            for r,v,w in [(src,aa,ns),(z,zz,L),(t,tt,nt)]:sim.write_across_shots(r,bits(v,w))
            sim.do(c)
            for r,v,w in [(src,aa,ns),(z,zz,L),(t,expected,nt)]:assert np.array_equal(sim.read_across_shots(r),bits(v,w)),(ns,L,nt,ex,sub)
            assert not sim.read_phase_flipped_across_shots().any()
            assert not sim.read_across_shots(km.array([km.q(i) for i in range(ns+L+nt,c.num_qubits)])).any()
            print('PASS Radix8',ns,L,nt,ex,sub,c.num_qubits,c.operation_statistics().num_toffoli_like,flush=True)
if __name__=='__main__':test()
