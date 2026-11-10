"""Complete baseline for the independently certified fused-square model.

This first version retains all polynomial coefficients through arithmetic.
It measures complete resources instead of assuming coefficient replay or
workspace scheduling is free. Host verifiers and the numeric model are unchanged.
"""
from pathlib import Path
import argparse,hashlib,json
import kickmix as km
import numpy as np
import affine_invsqrt as affine
import binary_centered_invsqrt as binary
import rounded_square
from booth_adder import accumulate
from short_adder import add
from monomial_unlookup import erase
from optimize_workspace import compact
from invsqrt import nr_reference_circuit as nr

HERE=Path(__file__).resolve().parent


def rounded_coordinate_square(b,x,h,capacity=32):
    b.x(x[13])
    for q in x[:13]:b.cx(x[13],q)
    guard=b.alloc_qubits(1);u=km.array([*x[5:13],guard[0]])
    add(b,km.array([x[4]]),u)
    rounded_square.apply(b,u,h,capacity)
    add(b,km.array([x[4]]),u,subtract=True);b.free(guard)
    for q in reversed(x[:13]):b.cx(x[13],q)
    b.x(x[13])


def lookup(b,address,registers,table):
    def load(row,control):
        for value,reg in zip(row,registers):b.gf2_iadd(int(value),target=reg,control=control)
    def visit(control,level,lo,hi):
        if np.all(table[lo:hi]==table[lo]):load(table[lo],control);return
        bit=address[level];mid=(lo+hi)//2;temporary=b.alloc_qubits(1)
        b.x(bit);b.init_and(control,bit,temporary[0]);b.x(bit)
        visit(temporary[0],level-1,lo,mid);b.cx(control,temporary[0]);visit(temporary[0],level-1,mid,hi)
        b.del_and(control,bit,temporary[0]);b.free(temporary)
    mid=len(table)//2;b.x(address[-1]);visit(address[-1],len(address)-2,0,mid);b.x(address[-1])
    visit(address[-1],len(address)-2,mid,len(table))


def restore_affine(b,x,y,flags,limit,inverse=False):
    head=y[28:]
    def constant(inv):
        a,c,d=y[24:27];temporary=head[0]
        if not inv:b.x(a);b.cx(a,c)
        b.init_and(a,c,temporary);b.x(temporary);b.cx(temporary,d);b.x(temporary);b.del_and(a,c,temporary)
        if inv:b.cx(a,c);b.x(a)
    def parity(inv):b.cx(flags[0],y[26])
    def coarse(inv):add(b,x[14:21],y[19:27],subtract=not inv,carry_workspace=head,qubit_limit=limit)
    operations=[constant,parity,coarse]
    for operation in (reversed(operations) if inverse else operations):operation(inverse)


def shift(b,y,flags,inverse=False):
    occupancy=[28]
    for k in range(4):occupancy.append(occupancy[-1]+(1<<k))
    for k in (reversed(range(4)) if inverse else range(4)):
        distance=1<<k;occ=occupancy[k];indices=[i for i in range(min(occ,len(y))) if i+distance<len(y)]
        for i in (indices if inverse else reversed(indices)):
            a,c=y[i],y[i+distance];control=flags[k+1]
            if i==27 and distance==1 and occ==28:b.cx(control,a);b.cx(control,c)
            elif i+distance>=occ:
                if inverse:b.cx(c,a);b.del_and(control,a,c)
                else:b.init_and(control,a,c);b.cx(c,a)
            else:b.cswap(control,a,c)


def build(limit=120,capacity=32):
    b=km.CircuitBuilder();x=b.create_quantum_register(22,'x');flags=b.alloc_qubits(5)
    def normalize(inverse=False):
        for k in (range(5) if inverse else reversed(range(5))):
            distance=1<<k
            if inverse:
                binary.cshift(b,x,distance,flags[k],22-distance,True)
                binary.del_flag_if_zero(b,x[-distance:],flags[k])
            else:
                b.flip_if_equal(x[-distance:],0,target=flags[k])
                binary.cshift(b,x,distance,flags[k],22-distance)
    normalize();h=b.alloc_qubits(11);rounded_coordinate_square(b,x,h,capacity)
    low=b.alloc_qubits(24);slope=b.alloc_qubits(17);curve=b.alloc_qubits(7)
    table=affine.table();table[:,0]<<=1;address=km.array([*x[13:21],flags[0]])
    lookup(b,address,[low,slope,curve],table)
    padding=b.alloc_qubits(13)
    allocated=km.array([*low,*padding]);y=km.array([*allocated[:27],x[21],*allocated[27:]])
    def curvature(inverse=False):
        b.x(x[21]);coordinate=km.array([*h,x[21]])
        accumulate(b,curve,coordinate,y[:27],-6,subtract=inverse,borrow=y[28:],qubit_limit=limit)
        b.x(x[21])
    def linear(inverse=False):
        b.x(x[13]);b.x(x[21])
        accumulate(b,slope,x[:14],y[:27],-11,subtract=not inverse,borrow=y[28:],auxiliary=km.array([x[21]]),qubit_limit=limit)
        b.x(x[21]);b.x(x[13])
    curvature();linear();restore_affine(b,x,y,flags,limit)
    shift(b,y,flags);nr._phase_by_bits(b,y);shift(b,y,flags,True)
    restore_affine(b,x,y,flags,limit,True);linear(True);curvature(True)
    packed=[int(row[0])|(int(row[1])<<24)|(int(row[2])<<41) for row in table]
    erase(b,address,km.array([*low,*slope,*curve]),packed)
    b.free(allocated);b.free(slope);b.free(curve)
    rounded_coordinate_square(b,x,h,capacity);b.free(h)
    normalize(True);b.free(flags);return b.finish_circuit()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=120);p.add_argument('--capacity',type=int,default=32);args=p.parse_args()
    path=HERE/'candidates'/f'fused_square_q{args.limit}_c{args.capacity}.kmb';c=compact(build(args.limit,args.capacity),path);c.to_file(str(path),format='kmb')
    row=dict(status='BUILT_REQUIRES_NATIVE_MODEL_CHECK',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    print(json.dumps(row),flush=True);(HERE/'results/fused_square_build.json').write_text(json.dumps(row,indent=2)+'\n')
