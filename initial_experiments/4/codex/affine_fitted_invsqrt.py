"""Screen newly fitted curvature widths using the selected affine schedule.

Production affine_invsqrt.build is copied with one guarded width substitution.
Every selected variant still requires the full native verifier.
"""
import argparse,inspect,json
from pathlib import Path
import numpy as np
import affine_invsqrt as affine
import binary_centered_invsqrt as p
from short_adder import add
from optimize_workspace import compact
import kickmix as km


def make_builder(raw,widths,correction):
    lookup=7 if np.array_equal(raw[::2],raw[1::2]) else 8
    tab=np.array(raw[::2] if lookup==7 else raw,dtype=np.int64);index=np.arange(len(tab));parity=index>>lookup;coarse=(index&((1<<lookup)-1))>>(lookup-7)
    tab[:,0]+=(1<<25)+(coarse<<18)+(1<<17)-np.array(affine.A)[parity]-correction
    assert tab.min()>=0 and int(tab[:,0].max()).bit_length()==23
    def restore(b,x,y,regs,flags,limit,inverse=False):
        def extra(inv):
            if not correction:return
            assert correction&(correction-1)==0
            bit=correction.bit_length()-1;head=y[27:];one=head[-1];b.x(one)
            add(b,km.array([one]),y[bit:27],subtract=inv,carry_workspace=head[:-1],qubit_limit=limit);b.x(one)
        if inverse:extra(True)
        affine.restore_affine(b,x,y,regs,flags,limit,inverse)
        if not inverse:extra(False)
    source=inspect.getsource(affine.build);old="namespace['build'](7,(27,18,8),";assert source.count(old)==1
    source=source.replace(old,"namespace['build'](7,FIT_WIDTHS,")
    assert source.count('lookup_frac=8')==1 and source.count('accum_bits=[18,9,None]')==1
    source=source.replace('lookup_frac=8','lookup_frac=FIT_LOOKUP').replace('accum_bits=[18,9,None]','accum_bits=FIT_SPLITS')
    namespace=dict(affine.__dict__);namespace.update(table=lambda:tab,FIT_WIDTHS=tuple(widths),FIT_LOOKUP=lookup,FIT_SPLITS=[18,widths[1]-9,None],restore_affine=restore)
    exec(source,namespace);return namespace['build']


def check_range(raw,widths):
    original=p.table_for;p.table_for=lambda *args:raw
    lo=2**62;hi=-2**62;worst=0.
    try:
        for begin in range(1,3133453,262144):
            xs=np.arange(begin,min(3133453,begin+262144),dtype=np.int64);k=21-np.floor(np.log2(xs)).astype(np.int64);s=k//2;normalized=xs<<k
            y=p.model(xs,7,widths,rounded=True,lookup_frac=8,booth=True)>>s
            coarse=(normalized>>14)&127;parity=k&1;residual=y+(1<<25)+(coarse<<18)+(1<<17)-np.array(affine.A)[parity]
            lo=min(lo,int(residual.min()));hi=max(hi,int(residual.max()))
            decoded=(y.astype(float)*np.exp2(s))/2**37;target=1/np.sqrt(xs);worst=max(worst,float(np.max(abs(decoded/target-1))))
    finally:p.table_for=original
    correction=0 if hi<1<<23 else 1<<17
    assert lo>=correction and hi-correction<1<<23,(lo,hi,correction)
    return dict(minimum=lo,maximum=hi,correction=correction,max_relative_error=worst)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--table',required=True);parser.add_argument('--w1',type=int,default=18);parser.add_argument('--w2',type=int,default=7);parser.add_argument('--tag',default='curvature7');parser.add_argument('--quick',action='store_true');a=parser.parse_args()
    raw=np.loadtxt(a.table,skiprows=1,dtype=np.int64);widths=(27,a.w1,a.w2);bounds=check_range(raw,widths);assert bounds['max_relative_error']<2**-24,bounds
    build=make_builder(raw,widths,bounds['correction']);rows=[]
    for limit in ([88,90,92,94] if a.quick else [86,88,89,90,92,94,96]):
        for keep in ([3,4] if a.quick else [2,3,4]):
            for delay in ([2] if a.quick else [2,3]):
                for product in ([17] if a.quick else [16,17]):
                    try:c=compact(build(limit=limit,keep=keep,delay=delay,product=product,deferred=3),Path(f'candidates/{a.tag}_screen_tmp.kmb'))
                    except (ValueError,AssertionError):continue
                    T=c.operation_statistics().num_toffoli_like;estimate=c.num_qubits*(T-11+251)
                    record=dict(Q=c.num_qubits,T=T,estimate=estimate,limit=limit,keep=keep,delay=delay,product=product)
                    if not rows or estimate<min(r['estimate'] for r in rows):c.to_file(f'candidates/{a.tag}_best.kmb',format='kmb')
                    rows.append(record)
    rows.sort(key=lambda r:r['estimate']);out=dict(widths=widths,bounds=bounds,rows=rows)
    Path(f'results/{a.tag}_resources.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(bounds=bounds,best=rows[:3])),flush=True)


if __name__=='__main__':main()
