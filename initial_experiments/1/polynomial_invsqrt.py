"""Direct piecewise minimax interpolation, with variable-precision Horner arithmetic."""
import argparse,json
from functools import lru_cache
from pathlib import Path
import numpy as np
import kickmix as km
from invsqrt import nr_reference_circuit as nr
from linear_invsqrt import qrom

@lru_cache(None)
def table_for(frac,widths,bias):
    order=len(widths)-1;h=2.**-frac
    table=np.zeros((2**(frac+2),order+1),dtype=np.int64)
    for i in range(2**frac,2**(frac+2)):
        a=i*h
        cheb=np.polynomial.Chebyshev.interpolate(lambda x:1/np.sqrt(a+x),order,domain=[0,h])
        coef=cheb.convert(kind=np.polynomial.Polynomial).coef
        for j,w in enumerate(widths):table[i,j]=round(coef[j]*(-1)**j*2**w)
    table[2**frac:,0]+=bias
    assert (table>=0).all()
    return table


def model(xs,frac,widths,bias=0):
    widths=tuple(widths);xs=np.asarray(xs,dtype=np.int64)
    s=(22-(np.floor(np.log2(xs)).astype(np.int64)+1))//2
    x=xs<<(2*s)
    regs=table_for(frac,widths,bias)[x>>(20-frac)].T.copy()
    for j in reversed(range(len(widths)-1)):
        for m in range(frac+1,21):
            e=widths[j]-widths[j+1]-m
            regs[j] -= ((x>>(20-m))&1)*(regs[j+1]<<e if e>=0 else regs[j+1]>>-e)
    return regs[0]<<s


def build(frac,widths,bias=0, lookup_delete=True, borrow=True, binary=True, dynamic=True):
    widths=tuple(widths);b=km.CircuitBuilder()
    x=b.create_quantum_register(22,'x')
    flags=b.alloc_qubits(4 if binary else 10)
    def normalize(inv):
        if not binary:return nr._normalize(b,inv,x,flags)
        for k in (range(4) if inv else reversed(range(4))):
            shift=2*(1<<k)
            if inv:b.right_rotate(x,shift,control=flags[k])
            b.flip_if_equal(x[-shift:],0,target=flags[k])
            if not inv:b.left_rotate(x,shift,control=flags[k])
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
            for m in range(frac+1,21):
                e=widths[j]-widths[j+1]-m
                terms=[(regs[j+1][i],i+e) for i in range(sizes[j+1]) if 0<=i+e<sizes[j]]
                if terms:
                    if dynamic:add_row(inv,x[20-m],terms,regs[j])
                    else:nr._add_row(b,inv,x[20-m],terms,regs[j],scratch,subtract=True)
    def shift(inv):
        if binary:
            for k in (range(4) if inv else reversed(range(4))):
                if inv:b.right_rotate(y,1<<k,control=flags[k])
                else:b.left_rotate(y,1<<k,control=flags[k])
            return
        for k in (reversed(range(10)) if inv else range(10)):
            lane=y[k:k+sizes[0]+1]
            if inv:b.right_rotate(lane,1,control=flags[k])
            else:b.left_rotate(lane,1,control=flags[k])
    qrom(b,x[20-frac:],regs,tab)
    horner(False);shift(False)
    nr._phase_by_bits(b,y[:widths[0]+10])
    shift(True);horner(True)
    if lookup_delete:
        packed=[sum(int(v) << sum(sizes[:j]) for j,v in enumerate(row)) for row in tab]
        b.del_lookup(table=packed,address=x[20-frac:],target=km.array([q for reg in regs for q in reg]))
    else:qrom(b,x[20-frac:],regs,tab)
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
