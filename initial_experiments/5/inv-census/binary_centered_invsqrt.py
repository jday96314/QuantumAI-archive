"""Binary normalization with parity folded into the coefficient lookup.

The mantissa always lies in [1,2); normalization parity contributes either
one or 1/sqrt(2) to the coefficient table, reducing the number of intervals.
"""
import argparse,json
from functools import lru_cache
from pathlib import Path
import numpy as np
import kickmix as km
from invsqrt import nr_reference_circuit as nr

@lru_cache(None)
def float_table(frac,order):
    h=2.**-frac
    table=np.zeros((2**(frac+1),order+1),dtype=float)
    for i in range(len(table)):
        a=1+(i%2**frac+0.5)*h
        scale=1 if i>>frac else 1/np.sqrt(2)
        cheb=np.polynomial.Chebyshev.interpolate(lambda x:scale/np.sqrt(a+x),order,domain=[-h/2,h/2])
        coef=cheb.convert(kind=np.polynomial.Polynomial).coef
        for j in range(order+1):table[i,j]=coef[j]*(-1)**j
    return table


@lru_cache(None)
def table_for(frac,widths,bias):
    table=np.rint(float_table(frac,len(widths)-1)*np.exp2(widths)).astype(np.int64)
    table[:,0]+=bias
    assert (table>=0).all()
    return table


def model(xs,frac,widths,bias=0,rounded=False,extra_frac=0,lookup_frac=None,rounding_masks=None,booth=False):
    widths=tuple(widths);xs=np.asarray(xs,dtype=np.int64)
    k=21-np.floor(np.log2(xs)).astype(np.int64)
    s=k//2;x=xs<<k
    lf=frac if lookup_frac is None else lookup_frac
    addr=((x>>(21-lf))&((1<<lf)-1))+((k&1)<<lf)
    regs=table_for(frac,widths,bias)[addr].T.copy()
    for j in reversed(range(len(widths)-1)):
        if booth is True or (isinstance(booth,(tuple,list)) and booth[j]):
            L=21-frac;z=(x&((1<<L)-1))^(1<<(L-1))
            for g in range((L+1)//2):
                bit=lambda k:0 if k<0 else (z>>min(k,L-1))&1
                digit=bit(2*g-1)+bit(2*g)-2*bit(2*g+1)
                for mag in [1,2]:
                    e=widths[j]-widths[j+1]-21+2*g+(mag==2)
                    value=regs[j+1]<<e if e>=0 else regs[j+1]>>-e
                    if rounded and e<0:value+=(regs[j+1]>>(-e-1))&1
                    if j==0:value=value<<extra_frac
                    regs[j]-=np.where(abs(digit)==mag,np.sign(digit)*value,0)
            continue
        for m in range(frac+1,22):
            e=widths[j]-widths[j+1]-m
            bit=(1-((x>>(21-m))&1)) if m==frac+1 else ((x>>(21-m))&1)
            value=regs[j+1]<<e if e>=0 else regs[j+1]>>-e
            if rounded and e<0 and (rounding_masks is None or ((rounding_masks[j]>>(m-frac-1))&1)):value+= (regs[j+1]>>(-e-1))&1
            if j==0:value=value<<extra_frac
            regs[j] += bit*value if m==frac+1 else -bit*value
    return regs[0]<<s


def cshift(b,reg,k,ctrl,occ,inv=False):
    """Controlled left shift by k of a value known to occupy reg[:occ] (inverse: right shift)."""
    idx=[i for i in range(min(occ,len(reg))) if i+k<len(reg)]
    for i in (idx if inv else reversed(idx)):b.cswap(ctrl,reg[i],reg[i+k])


def del_flag_if_zero(b,bits,flag):
    """Erases flag == (bits == 0) by X-basis measurement and a conditional phase fixup."""
    m=b.alloc_bits(1)[0]
    b.hmr(flag,m)
    b.push_condition(m)
    b.x(bits)
    if len(bits)==1:b.z(bits[0])
    elif len(bits)==2:b.cz(bits[0],bits[1])
    else:
        anc=b.alloc_qubits(len(bits)-2)
        b.init_and(bits[0],bits[1],anc[0])
        for i in range(1,len(bits)-2):b.init_and(anc[i-1],bits[i+1],anc[i])
        b.cz(anc[-1],bits[-1])
        for i in reversed(range(1,len(bits)-2)):b.del_and(anc[i-1],bits[i+1],anc[i])
        b.del_and(bits[0],bits[1],anc[0])
        b.free(anc)
    b.x(bits)
    b.pop_condition()


def build(frac,widths,bias=0, lookup_delete=True, borrow=True, binary=True, dynamic=True, accum_bits=None, carry_bits=2,rounded=False,extra_frac=0,global_phase=None,lookup_frac=None,reuse_normalized_top=False,output_extra_bits=11,curve_keep_bits=None,rounding_masks=None,booth=False,adder_prefix_bits=0,defer_workspace=False,qubit_limit=None,product_limit=None,defer_curve_phase=False,borrow_integration=False,parity_phase=False,monomial_delete=False,borrow_phase_headroom=False,unsigned_carries=False,delay_high_bits=0,monotone_delayed=False):
    widths=tuple(widths);b=km.CircuitBuilder()
    x=b.create_quantum_register(22,'x')
    assert binary
    flags=b.alloc_qubits(5)
    def normalize(inv):
        if not binary:return nr._normalize(b,inv,x,flags)
        for k in (range(5) if inv else reversed(range(5))):
            shift=1<<k
            if inv:
                cshift(b,x,shift,flags[k],22-shift,inv=True)
                del_flag_if_zero(b,x[-shift:],flags[k])
            else:
                b.flip_if_equal(x[-shift:],0,target=flags[k])
                cshift(b,x,shift,flags[k],22-shift)
    normalize(False)
    sizes=[max(1,int(tabcol.max()).bit_length()) for tabcol in table_for(frac,widths,bias).T]
    total_y=widths[0]+extra_frac+output_extra_bits-int(reuse_normalized_top)
    core_y=sizes[0]-int(reuse_normalized_top)
    if defer_workspace:assert dynamic
    allocated_y=b.alloc_qubits(core_y if defer_workspace else total_y)
    if reuse_normalized_top:
        assert np.all((table_for(frac,widths,bias)[:,0]>>(sizes[0]-1))==1)
        # Normalization leaves this input bit equal to one. It can carry the
        # constant top coefficient bit until Horner and shifting are reversed.
        y=km.array([*allocated_y[:sizes[0]-1],x[21],*allocated_y[sizes[0]-1:]])
    else:y=allocated_y
    regs=[y[:sizes[0]]]+[b.alloc_qubits(w) for w in sizes[1:]]
    work_regs=list(regs);work_regs[0]=regs[0][extra_frac:]
    arith_sizes=list(sizes);arith_sizes[0]-=extra_frac
    split=([accum_bits]+[None]*(len(regs)-1) if isinstance(accum_bits,int) else
           list(accum_bits) if accum_bits is not None else [None]*len(regs))
    assert len(split)==len(regs)
    counters=[] if defer_workspace else [b.alloc_qubits(carry_bits if low is not None else 0) for low in split]
    arith_regs=list(work_regs)
    for j,low in enumerate(split):
        if low is not None:
            assert 0<low<arith_sizes[j]
            if not defer_workspace:arith_regs[j]=km.array([*work_regs[j][:low],*counters[j]])
    scratch_extra=b.alloc_qubits(0 if dynamic else max(sizes)-(len(y)-sizes[0]) if borrow else max(sizes))
    scratch=km.array([*y[sizes[0]:],*scratch_extra]) if borrow else scratch_extra
    tab=table_for(frac,widths,bias)
    load_regs=list(regs)
    load_tab=tab.copy()
    if reuse_normalized_top:
        load_regs[0]=regs[0][:-1]
        load_tab[:,0]&=(1<<(sizes[0]-1))-1
    delayed=km.array([]);delayed_meas=None
    if delay_high_bits:
        assert reuse_normalized_top and monomial_delete and not borrow_phase_headroom and not unsigned_carries
        cutoff=len(load_regs[0])-delay_high_bits
        assert cutoff>=split[0]
        delayed=load_regs[0][cutoff:];load_regs[0]=load_regs[0][:cutoff]
        load_tab[:,0]&=(1<<cutoff)-1
    def workspace():return km.array([*y[sizes[0]:],*delayed])
    def add_row(inv,control,terms,target,round_source=None):
        low=min(i for q,i in terms);high=max(i for q,i in terms)
        assert [i for q,i in terms]==list(range(low,high+1))
        count=high-low+1;total=count+int(round_source is not None)
        borrowed=y[sizes[0]:sizes[0]+min(total,len(y)-sizes[0])] if borrow else km.array([])
        extra=b.alloc_qubits(total-len(borrowed))
        products=km.array([*borrowed,*extra]);sources=km.array([q for q,i in terms]+([round_source] if round_source is not None else []))
        b.init_and(control,sources,products)
        from short_adder import add
        add(b,products[:count],target[low:],subtract=not inv,carry_in=products[-1] if round_source is not None else None)
        b.del_and(control,sources,products);b.free(extra)
    def horner(inv):
        b.x(x[20-frac])
        steps=range(len(widths)-1) if inv else reversed(range(len(widths)-1))
        for j in steps:
            if inv and j==1 and curve_keep_bits is not None:restore_curve()
            if inv:integrate(j,True)
            use_booth=booth is True or (isinstance(booth,(tuple,list)) and booth[j])
            if use_booth:
                from booth_adder import accumulate
                # The leading normalized input bit is one between integrations.
                b.x(x[21])
                prefix=adder_prefix_bits if isinstance(adder_prefix_bits,int) else adder_prefix_bits[j]
                accumulate(b,work_regs[j+1],x[:21-frac],arith_regs[j],widths[j]-widths[j+1]-21,subtract=not inv,rounded=rounded,borrow=workspace() if borrow else (),auxiliary=km.array([x[21]]),adder_prefix_bits=prefix,qubit_limit=qubit_limit,product_limit=product_limit)
                b.x(x[21])
            for m in ([] if use_booth else range(frac+1,22)):
                e=widths[j]-widths[j+1]-m
                terms=[(work_regs[j+1][i],i+e) for i in range(arith_sizes[j+1]) if 0<=i+e<arith_sizes[j]]
                rs=regs[j+1][-e-1] if rounded and e<0 and -e-1<sizes[j+1] and (rounding_masks is None or ((rounding_masks[j]>>(m-frac-1))&1)) else None
                if terms:
                    if dynamic:add_row(not inv if m==frac+1 else inv,x[21-m],terms,arith_regs[j],rs)
                    else:nr._add_row(b,inv,x[21-m],terms,arith_regs[j],scratch,subtract=m!=frac+1)
                elif rs is not None:
                    # A rounded contribution can survive after its integer part vanishes.
                    product=b.alloc_qubits(1);b.init_and(x[21-m],rs,product[0])
                    from short_adder import add
                    add(b,product,arith_regs[j],subtract=not (not inv if m==frac+1 else inv))
                    b.del_and(x[21-m],rs,product[0]);b.free(product)
            if not inv:
                if j==0 and delay_high_bits:load_delayed()
                integrate(j,False)
            if not inv and j==1 and curve_keep_bits is not None:release_curve()
        b.x(x[20-frac])
    def shift(inv):
        if globals().get("NEW_PHASE_MODE"):return
        if binary:
            occ=[sizes[0]]
            for k in range(4):occ.append(occ[-1]+(1<<k))
            for k in (reversed(range(4)) if inv else range(4)):
                cshift(b,y,1<<k,flags[k+1],occ[k],inv=inv)
            return
        for k in (reversed(range(10)) if inv else range(10)):
            lane=y[k:k+sizes[0]+1]
            if inv:b.right_rotate(lane,1,control=flags[k])
            else:b.left_rotate(lane,1,control=flags[k])
    lf=frac if lookup_frac is None else lookup_frac
    address=km.array([*x[21-lf:21],*(() if parity_phase else (flags[0],))])
    def load_delayed():
        values=(tab[:,0]>>cutoff)&((1<<delay_high_bits)-1);anc=y[sizes[0]:]
        local_anc=b.alloc_qubits(max(0,len(address)-1-len(anc)))
        anc=km.array([*anc,*local_anc])
        if monotone_delayed:
            from monotone_lookup import prepare
            b.x(x[20-frac]);prepare(b,address,delayed,values,anc,x[21]);b.x(x[20-frac]);return
        assert len(anc)>=len(address)-1
        def visit(ctrl,level,lo,hi):
            if np.all(values[lo:hi]==values[lo]):b.gf2_iadd(int(values[lo]),target=delayed,control=ctrl);return
            bit=address[level];mid=(lo+hi)//2;t=anc[level]
            b.x(bit);b.init_and(ctrl,bit,t);b.x(bit)
            visit(t,level-1,lo,mid);b.cx(ctrl,t);visit(t,level-1,mid,hi);b.del_and(ctrl,bit,t)
        b.x(x[20-frac]);mid=len(values)//2;b.x(address[-1]);visit(address[-1],len(address)-2,0,mid);b.x(address[-1]);visit(address[-1],len(address)-2,mid,len(values));b.x(x[20-frac]);b.free(local_anc)
    curve_meas=None
    def release_curve():
        nonlocal curve_meas
        assert len(regs)==3 and 0<=curve_keep_bits<sizes[2]
        high=regs[2][curve_keep_bits:]
        # These unused phase-output bits are known zero here. Including them
        # supplies the deletion routine's local workspace without new wires.
        deletion_target=km.array([*high,*y[sizes[0]:]])
        if defer_curve_phase:
            curve_meas=b.alloc_bits(len(high));b.hmr(high,curve_meas)
        else:b.del_lookup(table=[int(v)>>curve_keep_bits for v in tab[:,2]],address=address,target=deletion_target)
        b.free(high)
    def restore_curve():
        nonlocal curve_meas
        high=b.alloc_qubits(sizes[2]-curve_keep_bits)
        regs[2]=km.array([*regs[2][:curve_keep_bits],*high])
        work_regs[2]=regs[2];arith_regs[2]=regs[2];load_regs[2]=regs[2]
        values=tab[:,2]>>curve_keep_bits
        def visit(ctrl,level,lo,hi):
            if np.all(values[lo:hi]==values[lo]):
                b.gf2_iadd(int(values[lo]),target=high,control=ctrl);return
            bit=address[level];mid=(lo+hi)//2;tmp=b.alloc_qubits(1)
            b.x(bit);b.init_and(ctrl,bit,tmp[0]);b.x(bit)
            visit(tmp[0],level-1,lo,mid);b.cx(ctrl,tmp[0]);visit(tmp[0],level-1,mid,hi)
            b.del_and(ctrl,bit,tmp[0]);b.free(tmp)
        # Horner holds the centered sign bit inverted here. The replay must
        # use the original lookup address, including that bit when present.
        if defer_curve_phase:b.x(x[20-frac])
        mid=len(values)//2;b.x(address[-1]);visit(address[-1],len(address)-2,0,mid);b.x(address[-1]);visit(address[-1],len(address)-2,mid,len(values))
        if defer_curve_phase:
            b.x(x[20-frac])
            # The measurement's phase is repaired when the same classical
            # function is available again; no separate phase lookup is needed.
            for cbit,qbit in zip(curve_meas,high):
                b.push_condition(cbit);b.z(qbit);b.pop_condition()
            b.bit_store(curve_meas,False);b.free(curve_meas);curve_meas=None
    def lookup():
        if lookup_frac is not None:
            def load(row,ctrl):
                for j,r in enumerate(load_regs):b.gf2_iadd(int(row[j]),target=r,control=ctrl)
            def visit(ctrl,level,lo,hi):
                if np.all(load_tab[lo:hi]==load_tab[lo]):load(load_tab[lo],ctrl);return
                bit=address[level];mid=(lo+hi)//2;work=b.alloc_qubits(1)
                b.x(bit);b.init_and(ctrl,bit,work[0]);b.x(bit)
                visit(work[0],level-1,lo,mid);b.cx(ctrl,work[0]);visit(work[0],level-1,mid,hi)
                b.del_and(ctrl,bit,work[0]);b.free(work)
            mid=len(tab)//2;b.x(address[-1]);visit(address[-1],len(address)-2,0,mid);b.x(address[-1]);visit(address[-1],len(address)-2,mid,len(tab));return
        with b.unary_iteration(address) as cursor:
            for a in range(len(tab)):
                cursor.move_to(a)
                for j,r in enumerate(load_regs):b.gf2_iadd(int(load_tab[a,j]),target=r,control=cursor.match_qubit)
    def integrate(j,inv):
        nonlocal delayed_meas
        if split[j] is None:return
        # horner temporarily flips the centered coordinate's sign bit.
        from short_adder import add
        high=work_regs[j][split[j]:arith_sizes[j]];counter=counters[j]
        if unsigned_carries:
            from signed_counter import add as add_carry
            b.x(x[20-frac])
            add_carry(b,counter,high,x[20-frac],subtract=inv,qubit_limit=qubit_limit,workspace=y[sizes[0]:] if borrow_integration else ())
            b.x(x[20-frac]);return
        if carry_bits==1:
            # Complement-conjugated increment is decrement. The centered
            # coordinate sign is currently inverted by horner(), so restore
            # it around the conjugation. This avoids a sign-extension AND.
            b.x(x[20-frac]);b.cx(x[20-frac],high);b.x(x[20-frac])
            signed=counter
        else:signed=km.array([counter[0],*[counter[1]]*(len(high)-1)])
        add(b,signed,high,subtract=inv,qubit_limit=qubit_limit,carry_workspace=(workspace() if j else y[sizes[0]:]) if borrow_integration else ())
        if carry_bits==1:
            b.x(x[20-frac]);b.cx(x[20-frac],high);b.x(x[20-frac])
        if inv and j==0 and delay_high_bits:
            delayed_meas=b.alloc_bits(len(delayed));b.hmr(delayed,delayed_meas)
    lookup()
    if defer_workspace:
        extra_y=b.alloc_qubits(0 if borrow_phase_headroom else total_y-core_y)
        allocated_y=km.array([*allocated_y,*extra_y]);y=km.array([*y,*extra_y])
        counters=[b.alloc_qubits(carry_bits if low is not None else 0) for low in split]
        for j,low in enumerate(split):
            if low is not None:arith_regs[j]=km.array([*work_regs[j][:low],*counters[j]])
    horner(False)
    if borrow_phase_headroom:
        assert defer_workspace and not parity_phase and len(y)==sizes[0]
        headroom=widths[0]+extra_frac+10-len(y)
        assert headroom<=10
        y=km.array([*y,*x[:headroom]])
    shift(False)
    phase_y=y[:widths[0]+extra_frac+10]
    if borrow_phase_headroom:
        from phase_headroom import apply
        apply(b,phase_y,sizes[0],flags[1:])
    else:
        mode=globals().get("NEW_PHASE_MODE")
        if mode=="borrowed":
            from exponent_kernels import borrowed_phase
            borrowed_phase(b,y[:27],x,flags)
        elif mode=="circular":
            from exponent_kernels import circular_phase
            circular_phase(b,y[:27],flags)
        elif mode=="walsh":
            from exponent_kernels import walsh_phase
            from fractions import Fraction
            import json
            co={int(k):Fraction(v) for k,v in json.loads((Path(__file__).resolve().parent/'results/walsh_search.json').read_text())['walsh_coefficients'].items()}
            walsh_phase(b,y[:27],flags[1:],co)
        else:nr._phase_by_bits(b,phase_y)
    if parity_phase:
        from fractions import Fraction
        correction=Fraction.from_float(float(np.sqrt(.5)))-1
        scratch_phase=b.alloc_qubits(1)[0]
        b.x(flags[0])
        for j,q in enumerate(phase_y):
            b.init_and(flags[0],q,scratch_phase)
            b.z_pow(scratch_phase,correction*Fraction(2)**(j-len(phase_y)+1))
            b.del_and(flags[0],q,scratch_phase)
        b.x(flags[0]);b.free(km.array([scratch_phase]))
    if global_phase:
        one=b.alloc_qubits(1);b.x(one);b.z_pow(one,2*global_phase);b.x(one);b.free(one)
    shift(True)
    if borrow_phase_headroom:y=y[:sizes[0]]
    horner(True)
    if defer_workspace:
        for counter in counters:b.free(counter)
        counters=[];b.free(extra_y);allocated_y=allocated_y[:core_y];y=y[:sizes[0]]
    if lookup_delete:
        load_sizes=[len(r) for r in load_regs]
        packed=[sum(int(v) << sum(load_sizes[:j]) for j,v in enumerate(row)) for row in load_tab]
        if delay_high_bits:
            packed=[v|((((int(row[0])>>cutoff)&((1<<delay_high_bits)-1)))<<sum(load_sizes)) for v,row in zip(packed,tab)]
        target=km.array([q for reg in load_regs for q in reg])
        if monomial_delete:
            from monomial_unlookup import erase
            erase(b,address,target,packed,extra_outcomes=delayed_meas if delay_high_bits else ())
        else:b.del_lookup(table=packed,address=address,target=target)
    else:lookup()
    normalize(True)
    for r in (allocated_y,flags,*regs[1:],scratch_extra,*counters):b.free(r)
    return b.finish_circuit()


def main():
    p=argparse.ArgumentParser();p.add_argument('--frac',type=int,default=6)
    p.add_argument('--widths',type=int,nargs='+',default=[29,23,17])
    p.add_argument('--bias',type=int,default=3)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();kw=dict(frac=a.frac,widths=tuple(a.widths),bias=a.bias)
    c=build(**kw);c.to_file(a.output,format='kmb')
    xs=np.arange(1,3133453);error=(model(xs,**kw)/2.**(a.widths[0]+10)*np.sqrt(xs)-1)
    print(json.dumps(dict(parameters=kw,qubits=c.num_qubits,toffolis=c.operation_statistics().num_toffoli_like,max_relative_error=float(np.max(np.abs(error))),min_signed_error=float(error.min()),max_signed_error=float(error.max()),worst_x=int(np.abs(error).argmax()+1)),indent=2))
if __name__=='__main__':main()
