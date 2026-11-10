"""Table-assisted linear initial estimate with reversible Newton refinement."""
import argparse
import json
from pathlib import Path
import numpy as np
import kickmix as km
from invsqrt import nr_reference_circuit as nr


def coefficients(frac, wc):
    # Uniform intervals in M, with a secant shifted to center its curvature error.
    a = np.arange(2**(frac+2), dtype=float) / 2**frac
    a = np.maximum(a, 1)
    h = 2.0**-frac
    slope = (1/np.sqrt(a) - 1/np.sqrt(a+h))/h
    intercept = 1/np.sqrt(a)
    miderr = intercept-slope*h/2-1/np.sqrt(a+h/2)
    intercept -= miderr/2
    return np.rint(np.stack([intercept, slope], axis=1)*2**wc).astype(np.int64)


def model(xs, frac, wc, wsq, wz, wp, bias=0):
    xs = np.asarray(xs, dtype=np.int64)
    s = (22-(np.floor(np.log2(xs)).astype(np.int64)+1))//2
    x = xs << (2*s)
    table = coefficients(frac,wc)
    c0,c1 = table[x >> (20-frac)].T
    y0 = c0.copy()
    for m in range(frac+1,21):
        if m < wc:
            y0 -= ((x >> (20-m)) & 1)*(c1 >> m)
    square = np.zeros_like(xs)
    for k in range(1,wc+1):
        if 2*k<=wsq: square += ((y0>>(wc-k))&1) << (wsq-2*k)
        for l in range(k+1,wc+1):
            if k+l-1<=wsq:
                square += (((y0>>(wc-k))&(y0>>(wc-l)))&1) << (wsq-k-l+1)
    z=np.zeros_like(xs)
    for i in range(22):
        e=wz-wsq+i-20
        z+=((x>>i)&1)*(square<<e if e>=0 else square>>-e)
    z=(3<<wz)-z
    p=np.zeros_like(xs)+bias
    for k in range(1,wc+1):
        e=wp-wz-k
        p+=((y0>>(wc-k))&1)*(z<<e if e>=0 else z>>-e)
    return p << s


def qrom(builder, address, regs, table):
    with builder.unary_iteration(address) as cursor:
        for a in range(1 << len(address)):
            if a < len(table) and a >= len(table)//4:
                cursor.move_to(a)
                for j,r in enumerate(regs):
                    builder.gf2_iadd(int(table[a,j]),target=r,control=cursor.match_qubit)


def build(frac=5,wc=17,wsq=28,wz=28,wp=29,bias=0):
    b=km.CircuitBuilder()
    x=b.create_quantum_register(22,'x')
    y1=b.alloc_qubits(wp+11)
    flags=b.alloc_qubits(10)
    coefs=[b.alloc_qubits(wc) for _ in range(2)]
    square=b.alloc_qubits(wsq)
    z=b.alloc_qubits(wz+2)
    scratch=b.alloc_qubits(max(wsq,wz+2,wp+1))
    table=coefficients(frac,wc)
    def initial(inv):
        def horner():
            for m in range(frac+1,min(20,wc-1)+1):
                nr._add_row(b, inv, x[20-m],[(coefs[1][i],i-m) for i in range(m,wc)],coefs[0],scratch,subtract=True)
        if inv: horner()
        qrom(b,x[20-frac:],coefs,table)
        if not inv: horner()
    def compute(inv):
        p=y1[:wp+1]
        if not inv:
            nr._normalize(b,False,x,flags)
            initial(False)
            nr._square(b,False,coefs[0],square,scratch)
            nr._multiply_m(b,False,x,square,z,scratch)
            nr._correct(b,False,z)
            nr._multiply_y0(b,False,coefs[0],z,p,scratch)
            if bias: b.iadd(bias,target=p)
            nr._denormalize(b,False,flags,y1)
        else:
            nr._denormalize(b,True,flags,y1)
            if bias: b.isub(bias,target=p)
            nr._multiply_y0(b,True,coefs[0],z,p,scratch)
            nr._correct(b,True,z)
            nr._multiply_m(b,True,x,square,z,scratch)
            nr._square(b,True,coefs[0],square,scratch)
            initial(True)
            nr._normalize(b,True,x,flags)
    compute(False)
    nr._phase_by_bits(b,y1)
    compute(True)
    for r in (y1,flags,*coefs,square,z,scratch): b.free(r)
    return b.finish_circuit()


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--frac',type=int,default=5)
    p.add_argument('--wc',type=int,default=17)
    p.add_argument('--wsq',type=int,default=28)
    p.add_argument('--wz',type=int,default=28)
    p.add_argument('--wp',type=int,default=29)
    p.add_argument('--bias',type=int,default=0)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    kw={k:v for k,v in vars(args).items() if k!='output'}
    circuit=build(**kw)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    circuit.to_file(args.output,format='kmb')
    xs=np.arange(1,3133453)
    e=np.abs(model(xs,**kw)/2.**(args.wp+11)*np.sqrt(xs)-1)
    print(json.dumps(dict(parameters=kw,qubits=circuit.num_qubits,toffolis=circuit.operation_statistics().num_toffoli_like,max_relative_error=float(e.max()),worst_x=int(e.argmax()+1)),indent=2))

if __name__=='__main__': main()
