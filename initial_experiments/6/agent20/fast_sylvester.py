"""Local source-division optimization for the read-only Euclidean Sylvester map.

Degree-one quotient division combines the two tail MACs into one bilinear
convolution. It preserves the prior armed-pivot convention and payload map.
Apache-2.0 arithmetic and prior algorithm are imported read-only.
"""
import importlib.util,sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
PRIOR=HERE/'../../5/dqi-core-phase'
if str(PRIOR.resolve()) not in sys.path:sys.path.append(str(PRIOR.resolve()))

def module(name,file):
    spec=importlib.util.spec_from_file_location(name,PRIOR/file)
    out=importlib.util.module_from_spec(spec);spec.loader.exec_module(out);return out

syl=module('agent20_original_sylvester','sylvester_circuit.py')
original_divide=syl.divide

def divide(m,E,f,g,inverse=False,one=False):
    n=len(g)-1;mm=len(f)-1
    if mm!=n+1 or n<3:return original_divide(m,E,f,g,inverse,one)
    q1,q0=f[-1],f[-2];tail=f[:-2]
    def shear():E.xor(m.prod_rt(q0,q1,g[-2]))
    def product():E.xor(m.polyprod_rt([q0,q1],g[:-1],[[w] for w in tail]+[[]]))
    with syl.factor(m,E,g[-1],one) as gi,m.opened(E,'fast_sylv_division_product',12) as tmp:
        def scale(w,back):
            if not one:m.mon_slot(E,w,g[-1] if back else gi,gi if back else g[-1],tmp)
        if inverse:product();scale(q0,True);shear();scale(q1,True)
        else:scale(q1,False);shear();scale(q0,False);product()

syl.divide=divide
apply=syl.apply
output_layout=syl.output_layout

