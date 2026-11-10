"""Replay curvature from the live coefficient rather than another table."""
import inspect,json,argparse
from pathlib import Path
import kickmix as km
import affine_invsqrt as affine
import binary_centered_invsqrt as poly
import fixed_top_affine_invsqrt as fixed
import zero_head_affine_invsqrt as zero
from optimize_workspace import compact


def change_source(source):
 source=source.replace('from booth_adder import accumulate','from booth_pair_carry_adder import accumulate')
 old='        mid=len(values)//2;b.x(address[-1]);visit(address[-1],len(address)-2,0,mid);b.x(address[-1]);visit(address[-1],len(address)-2,mid,len(values))'
 assert source.count(old)==1
 source=source.replace(old,"        from curve_register_restore import apply as replay_curve\n        replay_curve(b,regs[0],address,regs[2][:curve_keep_bits],high,curve_keep_bits)")
 source=source.replace('nonlocal curve_meas,allocated_y,y','nonlocal curve_meas,allocated_y,y,delayed_meas')
 old='        if inv and j==0 and delay_high_bits:\n            delayed_meas=b.alloc_bits(len(delayed));b.hmr(delayed,delayed_meas)'
 assert source.count(old)==1;source=source.replace(old,'')
 old='            b.bit_store(curve_meas,False);b.free(curve_meas);curve_meas=None'
 source=source.replace(old,old+'\n        if delay_high_bits:\n            delayed_meas=b.alloc_bits(len(delayed));b.hmr(delayed,delayed_meas)')
 return source
s=inspect.getsource(affine.build).replace("    (ROOT/'parameters/affine_build_source.py').write_text(source)",'    pass').replace('    namespace=dict(p.__dict__)','    source=change_source(source)\n    namespace=dict(p.__dict__)')
ns=dict(vars(affine));ns.update(restore_affine=fixed.restore_affine,change_source=change_source);exec(s,ns);poly.cshift=zero.cshift;build=ns['build']
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=87);p.add_argument('--product',type=int,default=17);p.add_argument('--keep',type=int,default=0);p.add_argument('--delay',type=int,default=2);p.add_argument('--deferred',type=int,default=5);a=p.parse_args();path=Path(f'candidates/pair_curve_q{a.limit}_p{a.product}_k{a.keep}_d{a.delay}_f{a.deferred}.kmb');c=compact(build(limit=a.limit,keep=a.keep,delay=a.delay,product=a.product,deferred=a.deferred),path);c.to_file(str(path),format='kmb');r=dict(status='BUILT',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,path=str(path));print(json.dumps(r),flush=True)
