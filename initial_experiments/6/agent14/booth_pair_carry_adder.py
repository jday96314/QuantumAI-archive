"""Radix-four multiplication with paired convolution products.

Two adjacent multiplier bits share one ripple addition. Partial products use
measurement erasure; the multiplier and multiplicand are restored exactly.
"""
import kickmix as km
from short_adder import add

def accumulate(b,source,coordinate,target,exponent,*,subtract=False,rounded=True,borrow=(),auxiliary=(),adder_prefix_bits=0,qubit_limit=None,product_limit=None):
    L=len(coordinate)
    for g in range((L+1)//2):
        top=2*g+1>=L
        e=exponent+2*g
        shapes=[]
        for branch in range(1 if top else 2):
            ee=e+branch
            shape=[(i,i+ee) for i in range(len(source)) if 0<=i+ee<len(target)]
            rs=-ee-1 if rounded and ee<0 and -ee-1<len(source) else None
            shapes.append((shape,rs))
        if not any(shape or rs is not None for shape,rs in shapes):continue
        current=coordinate[2*g];prev=coordinate[2*g-1] if g else None
        needed=1 if top else 2;extra_tmp=b.alloc_qubits(max(0,needed-len(auxiliary)))
        tmp=km.array([*auxiliary[:needed],*extra_tmp])
        if top:
            sign=tmp[0];b.cx(current,sign)
            if prev is not None:b.cx(prev,current)
            controls=[current]
        else:
            sign=coordinate[2*g+1];diff,two=tmp
            if prev is not None:b.cx(prev,current)
            b.cx(sign,diff);b.cx(current,diff)
            if prev is not None:b.cx(prev,diff)
            b.x(current);b.init_and(current,diff,two);b.x(current)
            controls=[current,two]
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
        if not top:
            if prev is not None:b.cx(prev,diff)
            b.cx(current,diff);b.cx(sign,diff)
        row_borrow=list(borrow)+([diff] if not top else [])
        borrowed=list(row_borrow[:min(total,len(row_borrow))]);extra=b.alloc_qubits(total-len(borrowed))
        products=km.array([*borrowed,*extra]);claims=[]
        for k in range(count):claims.append(normal.get(low+k,[]))
        if carry:claims.append(carry)
        pending=[list(terms) for terms in claims]
        def emit_product(c,s,targets):
            for out in targets[1:]:b.cx(targets[0],out)
            b.ccx(c,s,targets[0])
            for out in targets[1:]:b.cx(targets[0],out)
        if not top:
            # Each pair (s_i,s_{i+1}) times (c0+c1*x) needs three ANDs.
            # The middle coefficient combines A+B+C, with
            # C=(s_i+s_{i+1})*(c0+c1).
            where={low+k:k for k in range(count)}
            if carry:where[-1]=count
            i=0
            while i<len(source)-1:
                positions=[i+e,i+e+1,i+e+2]
                if not all(p in where for p in positions):i+=1;continue
                ids=[where[p] for p in positions];u,v=source[i:i+2];c,d=controls
                required=[(ids[0],(c,u)),(ids[1],(d,u)),(ids[1],(c,v)),(ids[2],(d,v))]
                if not all(term in pending[k] for k,term in required):i+=1;continue
                emit_product(c,u,[products[ids[0]],products[ids[1]]])
                emit_product(d,v,[products[ids[1]],products[ids[2]]])
                b.cx(d,c);b.cx(v,u);emit_product(c,u,[products[ids[1]]]);b.cx(v,u);b.cx(d,c)
                for k,term in required:pending[k].remove(term)
                i+=2
        for out,terms in zip(products,pending):
            for ctrl,src in terms:b.ccx(ctrl,src,out)
        lane=target[low:]
        b.cx(sign,lane)
        carry_candidates=([diff] if not top and diff not in products else [])+row_borrow[len(borrowed):]
        carry_work=km.array(list(dict.fromkeys(q for q in carry_candidates if q not in products)))
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
        if not top:
            b.cx(sign,diff);b.cx(current,diff)
            if prev is not None:b.cx(prev,diff)
        if top:
            if prev is not None:b.cx(prev,current)
            b.cx(current,sign)
        else:
            b.x(current);b.del_and(current,diff,two);b.x(current)
            if prev is not None:b.cx(prev,diff)
            b.cx(current,diff);b.cx(sign,diff)
            if prev is not None:b.cx(prev,current)
        b.free(extra_tmp)

def test():
    import numpy as np
    rng=np.random.default_rng(8102)
    for ns,L,nt,ex in [(5,5,8,-2),(5,6,9,-4),(4,5,7,0),(3,4,6,-6)]:
        for sub in [False,True]:
            b=km.CircuitBuilder();src=b.create_quantum_register(ns,'s');z=b.create_quantum_register(L,'z');t=b.create_quantum_register(nt,'t')
            accumulate(b,src,z,t,ex,subtract=sub);c=b.finish_circuit();N=4096
            aa=np.resize(np.repeat(np.arange(1<<ns),1<<L),N);zz=np.resize(np.tile(np.arange(1<<L),1<<ns),N);tt=rng.integers(0,1<<nt,N)
            delta=np.zeros(N,dtype=np.int64)
            for g in range((L+1)//2):
                bit=lambda k:0 if k<0 else (zz>>min(k,L-1))&1
                digit=bit(2*g-1)+bit(2*g)-2*bit(2*g+1)
                for mag in [1,2]:
                    e=ex+2*g+(mag==2);v=aa<<e if e>=0 else (aa>>-e)+((aa>>(-e-1))&1)
                    delta+=np.where(abs(digit)==mag,np.sign(digit)*v,0)
            expected=(tt-delta if sub else tt+delta)%(1<<nt)
            bits=lambda x,w:((x[None,:]>>np.arange(w)[:,None])&1).astype(bool)
            sim=km.Simulator(batch_size=N,seed=181);sim.use_same_registers_as(c)
            for r,v,w in [(src,aa,ns),(z,zz,L),(t,tt,nt)]:sim.write_across_shots(r,bits(v,w))
            sim.do(c)
            for r,v,w in [(src,aa,ns),(z,zz,L),(t,expected,nt)]:assert np.array_equal(sim.read_across_shots(r),bits(v,w)),(ns,L,nt,ex,sub)
            assert not sim.read_phase_flipped_across_shots().any()
            assert not sim.read_across_shots(km.array([km.q(i) for i in range(ns+L+nt,c.num_qubits)])).any()
            print('PASS Booth',ns,L,nt,ex,sub,c.num_qubits,c.operation_statistics().num_toffoli_like,flush=True)
if __name__=='__main__':test()
