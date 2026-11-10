"""Reversible components and measured live-state instrumentation."""
import collections
from fractions import Fraction
import numpy as np
import kickmix as km
from binary_centered_invsqrt import cshift,del_flag_if_zero
from short_adder import add

class Profile:
    def __init__(self):
        self.inner=km.CircuitBuilder();self.stage='start';self.peaks=collections.defaultdict(int);self.events=[]
    def __getattr__(self,name):return getattr(self.inner,name)
    def sample(self):self.peaks[self.stage]=max(self.peaks[self.stage],self.num_allocated_qubits)
    def alloc_qubits(self,n):
        r=self.inner.alloc_qubits(n);self.sample();return r
    def create_quantum_register(self,n,name):
        r=self.inner.create_quantum_register(n,name);self.sample();return r
    def mark(self,name):
        self.stage=name;self.sample();self.events.append(dict(stage=name,live=self.num_allocated_qubits));self.info_tag(name)

def normalize(b,x,flags,inverse=False):
    for k in (range(5) if inverse else reversed(range(5))):
        shift=1<<k
        if inverse:
            cshift(b,x,shift,flags[k],22-shift,inv=True);del_flag_if_zero(b,x[-shift:],flags[k])
        else:
            b.flip_if_equal(x[-shift:],0,target=flags[k]);cshift(b,x,shift,flags[k],22-shift)

def lookup(b,address,targets,rows):
    rows=np.asarray(rows,dtype=np.int64)
    if rows.ndim==1:rows=rows[:,None]
    def leaf(ctrl,row):
        for reg,value in zip(targets,row):b.gf2_iadd(int(value),target=reg,control=ctrl)
    def visit(ctrl,level,lo,hi):
        if np.all(rows[lo:hi]==rows[lo]):leaf(ctrl,rows[lo]);return
        assert level>=0
        bit=address[level];mid=(lo+hi)//2;t=b.alloc_qubits(1)[0]
        b.x(bit);b.init_and(ctrl,bit,t);b.x(bit)
        visit(t,level-1,lo,mid);b.cx(ctrl,t);visit(t,level-1,mid,hi)
        b.del_and(ctrl,bit,t);b.free(km.array([t]))
    mid=len(rows)//2;b.x(address[-1]);visit(address[-1],len(address)-2,0,mid)
    b.x(address[-1]);visit(address[-1],len(address)-2,mid,len(rows))

def shift(b,y,core,flags,inverse=False):
    """Zero-headroom shifts; reverses private moves by measured AND deletion."""
    stages=[];occ=core
    for j in range(4):stages.append((1<<j,occ,flags[j+1]));occ+=1<<j
    for step,occ,ctrl in (reversed(stages) if inverse else stages):
        idx=range(min(occ,len(y)-step))
        for i in (idx if inverse else reversed(idx)):
            src,dst=y[i],y[i+step]
            if i+step>=occ:
                if inverse:b.cx(dst,src);b.del_and(ctrl,src,dst)
                else:b.init_and(ctrl,src,dst);b.cx(dst,src)
            else:b.cswap(ctrl,src,dst)

def phase(b,y,denom=37,sign=None,coefficient=1):
    # sign conjugation: (v XOR all_ones*sign)-all_ones*sign=(-1)^sign*v.
    if sign is not None:b.cx(sign,y)
    for j,q in enumerate(y):b.z_pow(q,coefficient*Fraction(2)**(j-denom+1))
    if sign is not None:
        b.cx(sign,y);b.z_pow(sign,-coefficient*Fraction((1<<len(y))-1,1<<(denom-1)))

def integrate_slope(b,slope,counter,zsign,workspace,limit,inverse=False):
    high=slope[9:]
    b.x(zsign);b.cx(zsign,high);b.x(zsign)
    add(b,counter,high,subtract=inverse,carry_workspace=workspace,qubit_limit=limit)
    b.x(zsign);b.cx(zsign,high);b.x(zsign)
