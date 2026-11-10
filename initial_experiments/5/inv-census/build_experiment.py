"""Compile complete circuits with phase headroom absent during arithmetic."""
import argparse,json,math,time
from pathlib import Path
import kickmix as km
import binary_centered_invsqrt as p
import affine_invsqrt as a
from short_adder import add
from optimize_workspace import compact
H=Path(__file__).resolve().parent

def restore(b,x,y,regs,flags,limit,inverse=False):
    assert len(y)==27 and y[26]==x[21]
    scratch=b.alloc_qubits(7)
    def constant(inv):
        aa,cc,dd=y[23:26];t=scratch[0]
        if not inv:b.x(aa);b.cx(aa,cc)
        b.init_and(aa,cc,t);b.x(t);b.cx(t,dd);b.x(t);b.del_and(aa,cc,t)
        if inv:b.cx(aa,cc);b.x(aa)
    def parity(inv):b.cx(flags[0],y[25])
    def coarse(inv):add(b,x[14:21],y[18:26],subtract=not inv,carry_workspace=scratch,qubit_limit=limit)
    for f in (reversed([constant,parity,coarse]) if inverse else [constant,parity,coarse]):f(inverse)
    b.free(scratch)

def synth(r,tdirect):
    margin=19995185/2**64
    return (tdirect+r*(.56*math.log2(r/(math.pi*margin))+5.3))/4 if r else tdirect/4

def build(mode,limit,keep,delay,product):
    p.NEW_PHASE_MODE=mode;a.restore_affine=restore
    c=a.build(limit=limit,keep=keep,delay=delay,product=product,deferred=0)
    path=H/f'candidates/{mode}_q{limit}_k{keep}_d{delay}_p{product}.kmb'
    c=compact(c,path)
    st=c.operation_statistics();r=st.num_z_pow(min_clifford_hierarchy_level=4);td=st.num_t_like
    row=dict(mode=mode,limit=limit,keep=keep,delay=delay,product=product,Q=c.num_qubits,
        static_Toffolis=st.num_toffoli_like,rotations=r,direct_T=td,synthesis_work=synth(r,td),
        projected_mean_Toffolis=st.num_toffoli_like-11,
        projected_objective=c.num_qubits*(st.num_toffoli_like-11+synth(r,td)),path=f'candidates/{mode}_best.kmb')
    return c,row

def main():
    pa=argparse.ArgumentParser();pa.add_argument('--mode',choices=['borrowed','walsh','circular'],default='borrowed');pa.add_argument('--grid',action='store_true')
    args=pa.parse_args();best=math.inf;start=time.monotonic()
    grid=[(90,4,2,17)] if not args.grid else [(q,k,2,pr) for q in [72,76,80,84,88,90,94] for k in [0,2,4] for pr in [8,12,17]]
    with (H/f'results/{args.mode}_grid.jsonl').open('w') as out:
        for limit,keep,delay,product in grid:
            try:
                c,row=build(args.mode,limit,keep,delay,product)
                if row['projected_objective']<best:
                    best=row['projected_objective'];c.to_file(str(H/f'candidates/{args.mode}_best.kmb'),format='kmb')
                    (H/f'results/{args.mode}_best.json').write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(row),flush=True)
                row['status']='BUILT'
            except (ValueError,AssertionError) as ex:
                row=dict(status='INFEASIBLE',mode=args.mode,limit=limit,keep=keep,delay=delay,product=product,error=str(ex))
            out.write(json.dumps(row)+'\n');out.flush()
    print('seconds',time.monotonic()-start,flush=True)

if __name__=='__main__':main()
