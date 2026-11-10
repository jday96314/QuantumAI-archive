"""A coarser-address cubic Horner prototype with early cubic-word erasure."""
import inspect,argparse,json
from pathlib import Path
import numpy as np
import kickmix as km
import affine_invsqrt as affine
import binary_centered_invsqrt as poly
import fixed_top_affine_invsqrt as fixed
import zero_head_affine_invsqrt as zero
from short_adder import add
from optimize_workspace import compact
RAW=np.loadtxt('parameters/cubic_f6_19_10_5_fit.txt',dtype=np.int64,skiprows=1)


def table():
 T=RAW.copy();i=np.arange(len(T));T[:,0]+=(1<<25)+((i&127)>>1<<19)+(1<<17)-np.array(affine.A)[i>>7];assert T.min()>=0 and T[:,0].max()<1<<23;return T


def restore(b,x,y,regs,flags,limit,inverse=False):
 # Same modulo2^26 constant and parity restoration, with a six-bit coarse word.
 source=inspect.getsource(fixed.restore_affine).replace('x[14:21],y[18:26]','x[15:21],y[19:26]');ns=dict(vars(fixed));exec(source,ns);ns['restore_affine'](b,x,y,regs,flags,limit,inverse)


def build(limit=90,product=17,keep=4,deferred=3):
 s=Path('../../4/codex/parameters/affine_build_source.py').read_text().replace('from booth_adder import accumulate','from booth_pair_carry_adder import accumulate').replace('assert len(regs)==3','assert len(regs)==4')
 # Parameterize the two existing deferred headroom allocations consistently.
 s=s.replace('total_y-core_y-3','total_y-core_y-'+str(deferred)).replace('b.alloc_qubits(3)\n        allocated_y','b.alloc_qubits('+str(deferred)+')\n        allocated_y').replace('allocated_y=allocated_y[:-3];y=y[:-3]','allocated_y=allocated_y[:-'+str(deferred)+'];y=y[:-'+str(deferred)+']')
 old='''        extra_y=b.alloc_qubits(0 if borrow_phase_headroom else total_y-core_y-'''+str(deferred)+''')
        allocated_y=km.array([*allocated_y,*extra_y]);y=km.array([*allocated_y[:26],x[21],*allocated_y[26:]])'''
 assert s.count(old)==1;s=s.replace(old,'        extra_y=None')
 s=s.replace('b.free(extra_y);allocated_y=allocated_y[:core_y]','allocated_y=allocated_y[:core_y]')
 s=s.replace('            if inv and j==1 and curve_keep_bits is not None:restore_curve()','            if inv and j==2:restore_cubic()\n            if inv and j==1 and curve_keep_bits is not None:restore_curve()')
 s=s.replace('            if not inv and j==1 and curve_keep_bits is not None:release_curve()','            if not inv and j==2:release_cubic()\n            if not inv and j==1 and curve_keep_bits is not None:release_curve()')
 a=s.index('    def release_curve():')
 s=s[:a]+'''    cubic_meas=None
    def release_cubic():
        nonlocal cubic_meas,extra_y,allocated_y,y
        cubic_meas=b.alloc_bits(len(regs[3]));b.hmr(regs[3],cubic_meas);b.free(regs[3])
        extra_y=b.alloc_qubits(total_y-core_y-DEFERRED)
        allocated_y=km.array([*allocated_y,*extra_y]);y=km.array([*allocated_y[:26],x[21],*allocated_y[26:]])
    def restore_cubic():
        nonlocal extra_y,allocated_y,y,cubic_meas
        b.free(extra_y);allocated_y=allocated_y[:core_y];y=y[:sizes[0]]
        word=b.alloc_qubits(sizes[3]);regs[3]=word;work_regs[3]=word;arith_regs[3]=word;load_regs[3]=word
        values=tab[:,3]
        def visit(ctrl,level,lo,hi):
            if np.all(values[lo:hi]==values[lo]):b.gf2_iadd(int(values[lo]),target=word,control=ctrl);return
            bit=address[level];mid=(lo+hi)//2;q=b.alloc_qubits(1)[0];b.x(bit);b.init_and(ctrl,bit,q);b.x(bit);visit(q,level-1,lo,mid);b.cx(ctrl,q);visit(q,level-1,mid,hi);b.del_and(ctrl,bit,q);b.free(km.array([q]))
        b.x(x[20-frac]);mid=len(values)//2;b.x(address[-1]);visit(address[-1],len(address)-2,0,mid);b.x(address[-1]);visit(address[-1],len(address)-2,mid,len(values));b.x(x[20-frac])
        for bit,q in zip(cubic_meas,word):b.push_condition(bit);b.z(q);b.pop_condition()
        b.bit_store(cubic_meas,False);b.free(cubic_meas);cubic_meas=None
'''.replace('DEFERRED',str(deferred))+s[a:]
 ns=dict(vars(poly));ns.update(table_for=lambda *args:table(),restore_affine=restore,cshift=zero.cshift);exec(s,ns)
 return ns['build'](6,(27,19,10,5),lookup_frac=7,accum_bits=[18,10,None,None],carry_bits=1,rounded=True,booth=True,reuse_normalized_top=True,output_extra_bits=10,defer_workspace=True,qubit_limit=limit,curve_keep_bits=keep,defer_curve_phase=True,product_limit=product,borrow_integration=True,monomial_delete=True,delay_high_bits=2)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=92);p.add_argument('--product',type=int,default=17);p.add_argument('--keep',type=int,default=4);p.add_argument('--deferred',type=int,default=3);a=p.parse_args();out=Path(f'candidates/cubic_horner_q{a.limit}_p{a.product}_k{a.keep}_d{a.deferred}.kmb');c=compact(build(a.limit,a.product,a.keep,a.deferred),out);c.to_file(str(out),format='kmb');print(json.dumps(dict(Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,path=str(out))),flush=True)
