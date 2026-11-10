import json,numpy as np,kickmix as km
import polynomial_invsqrt as poly, linear_invsqrt, short_adder
from invsqrt import nr_reference_circuit as nr
params=json.load(open('parameters/solution.json'))['invsqrt']
table=np.load('parameters/invsqrt_table.npy');orig=poly.table_for
poly.table_for=lambda f,w,b=0: table
def count(patch={}):
    saved={}
    for (mod,name),fn in patch.items(): saved[(mod,name)]=getattr(mod,name); setattr(mod,name,fn)
    try: c=poly.build(params['frac'],tuple(params['widths']),0)
    finally:
        for (mod,name),fn in saved.items(): setattr(mod,name,fn)
    return c.operation_statistics().num_toffoli_like, c.num_qubits
full=count()
print('full',full)
print('no qrom',count({(poly,'qrom'):lambda *a,**k:None}))
print('no phase',count({(nr,'_phase_by_bits'):lambda *a,**k:None}))
print('no adder',count({(short_adder,'add'):lambda *a,**k:None}))
B=km.CircuitBuilder
class NoRot(B): 
    def right_rotate(self,*a,**k): pass
    def left_rotate(self,*a,**k): pass
print('no rotates',count({(km,'CircuitBuilder'):NoRot}))
class NoDel(B):
    def del_lookup(self,*a,**k): pass
print('no del_lookup',count({(km,'CircuitBuilder'):NoDel}))
class NoEq(B):
    def flip_if_equal(self,*a,**k): pass
print('no flip_if_equal',count({(km,'CircuitBuilder'):NoEq}))
