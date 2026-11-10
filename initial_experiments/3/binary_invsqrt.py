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
from linear_invsqrt import qrom

@lru_cache(None)
def float_table(frac,order):
    h=2.**-frac
    table=np.zeros((2**(frac+1),order+1),dtype=float)
    for i in range(len(table)):
        a=1+(i%2**frac)*h
        scale=1 if i>>frac else 1/np.sqrt(2)
        cheb=np.polynomial.Chebyshev.interpolate(lambda x:scale/np.sqrt(a+x),order,domain=[0,h])
        coef=cheb.convert(kind=np.polynomial.Polynomial).coef
        for j in range(order+1):table[i,j]=coef[j]*(-1)**j
    return table


@lru_cache(None)
def table_for(frac,widths,bias):
    table=np.rint(float_table(frac,len(widths)-1)*np.exp2(widths)).astype(np.int64)
    table[:,0]+=bias
    assert (table>=0).all()
    return table


def model(xs,frac,widths,bias=0):
    widths=tuple(widths);xs=np.asarray(xs,dtype=np.int64)
    k=21-np.floor(np.log2(xs)).astype(np.int64)
    s=k//2;x=xs<<k
    addr=((x>>(21-frac))&((1<<frac)-1))+((k&1)<<frac)
    regs=table_for(frac,widths,bias)[addr].T.copy()
    for j in reversed(range(len(widths)-1)):
        for m in range(frac+1,22):
            e=widths[j]-widths[j+1]-m
            bit=(x>>(21-m))&1
            value=regs[j+1]<<e if e>=0 else regs[j+1]>>-e
            regs[j] -= bit*value
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


def build(frac,widths,bias=0, lookup_delete=True, borrow=True, binary=True, dynamic=True):
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
    y=b.alloc_qubits(widths[0]+11)
    sizes=[max(1,int(tabcol.max()).bit_length()) for tabcol in table_for(frac,widths,bias).T]
    regs=[y[:sizes[0]]]+[b.alloc_qubits(w) for w in sizes[1:]]
    scratch_extra=b.alloc_qubits(0 if dynamic else max(sizes)-(len(y)-sizes[0]) if borrow else max(sizes))
    scratch=km.array([*y[sizes[0]:],*scratch_extra]) if borrow else scratch_extra
    tab=table_for(frac,widths,bias)
    def add_row(inv,control,terms,target):
        low=min(i for q,i in terms);high=max(i for q,i in terms)
        assert [i for q,i in terms]==list(range(low,high+1))
        count=high-low+1
        borrowed=y[sizes[0]:sizes[0]+min(count,len(y)-sizes[0])] if borrow else km.array([])
        extra=b.alloc_qubits(count-len(borrowed))
        products=km.array([*borrowed,*extra]);sources=km.array([q for q,i in terms])
        b.init_and(control,sources,products)
        from short_adder import add
        add(b,products,target[low:],subtract=not inv)
        b.del_and(control,sources,products);b.free(extra)
    def horner(inv):
        steps=range(len(widths)-1) if inv else reversed(range(len(widths)-1))
        for j in steps:
            for m in range(frac+1,22):
                e=widths[j]-widths[j+1]-m
                terms=[(regs[j+1][i],i+e) for i in range(sizes[j+1]) if 0<=i+e<sizes[j]]
                if terms:
                    if dynamic:add_row(inv,x[21-m],terms,regs[j])
                    else:nr._add_row(b,inv,x[21-m],terms,regs[j],scratch,subtract=True)
    def shift(inv):
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
    address=km.array([*x[21-frac:21],flags[0]])
    def lookup():
        with b.unary_iteration(address) as cursor:
            for a in range(len(tab)):
                cursor.move_to(a)
                for j,r in enumerate(regs):b.gf2_iadd(int(tab[a,j]),target=r,control=cursor.match_qubit)
    lookup()
    horner(False);shift(False)
    nr._phase_by_bits(b,y[:widths[0]+10])
    shift(True);horner(True)
    if lookup_delete:
        packed=[sum(int(v) << sum(sizes[:j]) for j,v in enumerate(row)) for row in tab]
        b.del_lookup(table=packed,address=address,target=km.array([q for reg in regs for q in reg]))
    else:lookup()
    normalize(True)
    for r in (y,flags,*regs[1:],scratch_extra):b.free(r)
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
