"""Research-only affine preconditioning of the existing passing polynomial.

Store c0 + center/4 - A(parity); keep the other coefficients unchanged.
Restore A(parity) - center/4 before phasing, using only seven address bits.
The source generator below makes checked
replacements in a copy of the production builder, leaving it unchanged.
"""
import argparse,inspect,json
from pathlib import Path
import numpy as np
import kickmix as km
import binary_centered_invsqrt as p
from optimize_workspace import compact
from short_adder import add

ROOT=Path(__file__).resolve().parent
A=(961<<17,1217<<17)


def table():
    t=np.load(ROOT/'candidates/margin_curve8_t5.9e-08.npy').copy()
    index=np.arange(len(t));parity=index>>8;coarse=(index&255)>>1
    t[:,0] += (1<<25)+(coarse<<18)+(1<<17)-np.array(A)[parity]
    assert int(t[:,0].max()).bit_length()==23 and t[:,0].min()>=0
    return t


def restore_affine(b,x,y,regs,flags,limit,inverse=False):
    head=y[27:]
    assert len(head)==10 and x[21]==y[26]
    def constant(inv):
        one=head[-1];b.x(one)
        positions=[23,24]
        for bit in (reversed(positions) if inv else positions):
            add(b,km.array([one]),y[bit:27],subtract=inv,carry_workspace=head[:-1],qubit_limit=limit)
        b.x(one)
    def parity(inv):
        add(b,km.array([flags[0]]),y[25:27],subtract=inv,carry_workspace=head,qubit_limit=limit)
    def coarse(inv):
        add(b,x[14:21],y[18:27],subtract=not inv,carry_workspace=head,qubit_limit=limit)
    steps=[constant,parity,coarse]
    for step in (reversed(steps) if inverse else steps):step(inverse)


def build(limit=88,keep=4,delay=4,product=17,deferred=3):
    source=inspect.getsource(p.build)
    def replace(old,new):
        nonlocal source
        assert source.count(old)==1,(old,source.count(old))
        source=source.replace(old,new)
    replace('core_y=sizes[0]-int(reuse_normalized_top)','core_y=sizes[0]')
    replace('if reuse_normalized_top:\n        assert np.all','if False:\n        assert np.all')
    replace('if reuse_normalized_top:\n        load_regs[0]','if False:\n        load_regs[0]')
    replace('    regs=[y[:sizes[0]]]','    def spare():return km.array([q for q in y[sizes[0]:] if q!=x[21]])\n    regs=[y[:sizes[0]]]')
    # Longer free headroom replaces the omitted c0 bits; the normalized top
    # input bit joins the phase word only after the coefficient lookup.
    replace('allocated_y=km.array([*allocated_y,*extra_y]);y=km.array([*y,*extra_y])',
            'allocated_y=km.array([*allocated_y,*extra_y]);y=km.array([*allocated_y[:26],x[21],*allocated_y[26:]])')
    replace('total_y-core_y)\n        allocated_y=', f'total_y-core_y-{deferred})\n        allocated_y=')
    replace('curve_meas=None\n    def release_curve():', 'curve_meas=None\n    affine_extra=None\n    def release_curve():')
    replace('    def release_curve():\n        nonlocal curve_meas',
            '    def release_curve():\n        nonlocal curve_meas,affine_extra,allocated_y,y')
    replace('        b.free(high)\n    def restore_curve():',
            f'        b.free(high)\n        affine_extra=b.alloc_qubits({deferred})\n        allocated_y=km.array([*allocated_y,*affine_extra]);y=km.array([*y,*affine_extra])\n    def restore_curve():')
    replace('    def restore_curve():\n        nonlocal curve_meas',
            f'    def restore_curve():\n        nonlocal curve_meas,allocated_y,y\n        b.free(affine_extra);allocated_y=allocated_y[:-{deferred}] if {deferred} else allocated_y;y=y[:-{deferred}] if {deferred} else y')
    replace('occ=[sizes[0]]','occ=[widths[0]]')
    replace('horner(False)\n    if borrow_phase_headroom:',
            'horner(False)\n    restore_affine(b,x,y,regs,flags,qubit_limit)\n    if borrow_phase_headroom:')
    replace('shift(True)\n    if borrow_phase_headroom:',
            'shift(True)\n    restore_affine(b,x,y,regs,flags,qubit_limit,True)\n    if borrow_phase_headroom:')
    source=source.replace('y[sizes[0]:sizes[0]+min(total,len(y)-sizes[0])]','spare()[:min(total,len(spare()))]')
    source=source.replace('y[sizes[0]:]','spare()')
    # Avoid replacing the implementation of spare() recursively.
    source=source.replace('def spare():return km.array([q for q in spare() if q!=x[21]])',
                          'def spare():return km.array([q for q in y[sizes[0]:] if q!=x[21]])')
    namespace=dict(p.__dict__);namespace.update(table_for=lambda *args:table(),restore_affine=restore_affine)
    exec(source,namespace)
    (ROOT/'parameters/affine_build_source.py').write_text(source)
    return namespace['build'](7,(27,18,8),lookup_frac=8,accum_bits=[18,9,None],carry_bits=1,
          rounded=True,booth=True,reuse_normalized_top=True,output_extra_bits=0 if getattr(p,"NEW_PHASE_MODE",None) else 10,defer_workspace=True,
          qubit_limit=limit,curve_keep_bits=keep,defer_curve_phase=True,product_limit=product,
          borrow_integration=True,monomial_delete=True,delay_high_bits=delay)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--limit',type=int,default=88)
    parser.add_argument('--keep',type=int,default=4);parser.add_argument('--delay',type=int,default=4)
    args=parser.parse_args();path=ROOT/f'candidates/affine_q{args.limit}_k{args.keep}_d{args.delay}.kmb'
    c=compact(build(args.limit,args.keep,args.delay),path);c.to_file(str(path),format='kmb')
    print(json.dumps(dict(Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,path=str(path))),flush=True)
