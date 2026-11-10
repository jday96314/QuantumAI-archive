"""Resource-screen guarded C0 precision with a narrow curvature table.

Research only: all native inputs must be checked before selecting a circuit.
The added precision shifts the affine remainder and its restoration by one bit.
"""
import argparse, inspect, importlib.util, json
from pathlib import Path
import numpy as np
import kickmix as km
import affine_invsqrt as affine
import binary_centered_invsqrt as poly
from affine8_fitted_invsqrt import padded_erase
from optimize_workspace import compact
from short_adder import add


def make_builder(raw,w1,w2):
    lookup=(len(raw)//2).bit_length()-1
    if lookup>7 and np.array_equal(raw[::2],raw[1::2]):raw=raw[::2];lookup-=1
    idx=np.arange(len(raw));parity=idx>>lookup;coarse=(idx&((1<<lookup)-1))>>(lookup-7)
    table=raw.copy();table[:,0]+=(coarse<<19)-(11<<24)-(parity<<26)
    assert table.min()>=0 and int(table[:,0].max()).bit_length()<=24
    def restore(b,x,y,regs,flags,limit,inverse=False):
        head=y[28:];assert len(head)==10 and x[21]==y[27]
        def constant(inv):
            one=head[-1];b.x(one)
            for bit in (list(reversed([24,25])) if inv else [24,25]):
                add(b,km.array([one]),y[bit:28],subtract=inv,carry_workspace=head[:-1],qubit_limit=limit)
            b.x(one)
        def parity_restore(inv):add(b,km.array([flags[0]]),y[26:28],subtract=inv,carry_workspace=head,qubit_limit=limit)
        def coarse_restore(inv):add(b,x[14:21],y[19:28],subtract=not inv,carry_workspace=head,qubit_limit=limit)
        for step in (list(reversed([constant,parity_restore,coarse_restore])) if inverse else [constant,parity_restore,coarse_restore]):step(inverse)
    src=inspect.getsource(affine.build)
    changes={"namespace['build'](7,(27,18,8),":"namespace['build'](7,FIT_WIDTHS,",
             'lookup_frac=8':'lookup_frac=FIT_LOOKUP',
             'accum_bits=[18,9,None]':'accum_bits=[19,9,None]',
             '*allocated_y[:26],x[21],*allocated_y[26:]':'*allocated_y[:27],x[21],*allocated_y[27:]'}
    for old,new in changes.items():assert src.count(old)==1,(old,src.count(old));src=src.replace(old,new)
    ps=inspect.getsource(poly.build)
    old='erase(b,address,target,packed,extra_outcomes=delayed_meas if delay_high_bits else ())'
    assert ps.count(old)==1;ps=ps.replace(old,'padded_erase(b,address,target,packed,extra_outcomes=delayed_meas if delay_high_bits else ())')
    path=Path(f'parameters/guard28_poly_source_{w1}_{w2}.py');path.write_text(ps)
    spec=importlib.util.spec_from_file_location('guard28_poly',path);mod=importlib.util.module_from_spec(spec)
    mod.__dict__.update(poly.__dict__);mod.padded_erase=padded_erase
    exec(compile(ps,str(path.resolve()),'exec'),mod.__dict__)
    ns=dict(affine.__dict__);ns.update(table=lambda:table,restore_affine=restore,p=mod,FIT_WIDTHS=(28,w1,w2),FIT_LOOKUP=lookup)
    exec(src,ns);return ns['build']


def main():
    p=argparse.ArgumentParser();p.add_argument('--table',required=True);p.add_argument('--w1',type=int,default=18);p.add_argument('--w2',type=int,default=6);p.add_argument('--tag',default='guard28_18_6');a=p.parse_args()
    raw=np.loadtxt(a.table,skiprows=1,dtype=np.int64);build=make_builder(raw,a.w1,a.w2);rows=[]
    for limit in [82,84,86,88,90,92]:
        for keep in [1,2,3,4]:
            for product in [17,None]:
                try:c=compact(build(limit=limit,keep=keep,delay=2,product=product,deferred=3),Path(f'candidates/{a.tag}_tmp.kmb'))
                except (AssertionError,ValueError):continue
                T=c.operation_statistics().num_toffoli_like;estimate=c.num_qubits*(T-11+251)
                row=dict(limit=limit,keep=keep,product=product,Q=c.num_qubits,T=T,estimate=estimate)
                if not rows or estimate<min(r['estimate'] for r in rows):c.to_file(f'candidates/{a.tag}_best.kmb',format='kmb');print(json.dumps(row),flush=True)
                rows.append(row)
    rows.sort(key=lambda r:r['estimate']);Path(f'results/{a.tag}_resources.json').write_text(json.dumps(dict(widths=[28,a.w1,a.w2],rows=rows),indent=2)+'\n')
    print(json.dumps(dict(status='FINISHED',best=rows[:3])),flush=True)
if __name__=='__main__':main()
