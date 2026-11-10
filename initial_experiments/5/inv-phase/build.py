"""Complete incremental phase circuits, with the incumbent table unchanged."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import kickmix as km
from booth_adder import accumulate
from streamed_unlookup import erase
from optimize_workspace import compact
from components import Profile,normalize,lookup,shift,phase,integrate_slope
from model import HERE,TABLE,synthesis

def build(mode='split',limit=84,keep=0):
    b=Profile();x=b.create_quantum_register(22,'x');flags=b.alloc_qubits(5)
    b.mark('normalize');normalize(b,x,flags)
    address=km.array([*x[13:21],flags[0]])
    # Bias lives in the phase table, not in an extra phase correction.
    constants=TABLE[:,0]-(1<<18) if mode=='split' else TABLE[:,0]
    b.mark('constant_lookup');constant=b.alloc_qubits(37)
    lookup(b,address,[constant[:27]],constants)
    b.mark('constant_phase');shift(b,constant,27,flags);phase(b,constant);shift(b,constant,27,flags,True)
    constant_meas=b.alloc_bits(27);b.hmr(constant[:27],constant_meas);b.free(constant)
    b.mark('derivative_lookup');slope=b.alloc_qubits(17);curve=b.alloc_qubits(7)
    lookup(b,address,[slope,curve],TABLE[:,1:])
    counter=b.alloc_qubits(1)
    # The coordinate stays in the input wires; its sign is inverted for z.
    b.x(x[13])
    workspace=b.alloc_qubits(10)
    def curve_horner(inverse=False):
        if inverse:integrate_slope(b,slope,counter,x[13],workspace,limit,True)
        b.x(x[21])
        accumulate(b,curve,x[:14],km.array([*slope[:9],*counter]),-11,
            subtract=not inverse,rounded=True,borrow=workspace,auxiliary=km.array([x[21]]),qubit_limit=limit)
        b.x(x[21])
        if not inverse:integrate_slope(b,slope,counter,x[13],workspace,limit)
    b.mark('curve_horner');curve_horner()
    b.free(workspace)
    high=curve[keep:];curve_meas=b.alloc_bits(len(high));b.hmr(high,curve_meas);b.free(high)
    b.mark('residual_phase')
    if mode=='split':
        y=b.alloc_qubits(29);b.x(y[18])
        def residual(inverse=False):
            b.x(x[21])
            accumulate(b,slope,x[:14],y[:19],-12,subtract=not inverse,
                rounded=True,borrow=y[19:],auxiliary=km.array([x[21]]),qubit_limit=limit)
            b.x(x[21])
        residual();shift(b,y,19,flags);phase(b,y);shift(b,y,19,flags,True);residual(True)
        b.x(y[18]);b.free(y)
    elif mode=='rows':
        from direct_rows import apply
        apply(b,slope,x[:14],flags,-12)
    else:raise ValueError(mode)
    b.mark('curvature_replay');new_high=b.alloc_qubits(7-keep);curve=km.array([*curve[:keep],*new_high])
    b.x(x[13]);lookup(b,address,[new_high],TABLE[:,2]>>keep);b.x(x[13])
    for bit,q in zip(curve_meas,new_high):
        b.push_condition(bit);b.z(q);b.pop_condition()
    b.bit_store(curve_meas,False);b.free(curve_meas)
    workspace=b.alloc_qubits(10);b.mark('reverse_curve_horner');curve_horner(True);b.free(workspace)
    b.free(counter);b.x(x[13])
    b.mark('table_phase_repair')
    packed=[int(row[1])|(int(row[2])<<17)|(int(c) <<24) for row,c in zip(TABLE,constants)]
    erase(b,address,km.array([*slope,*curve]),packed,extra_outcomes=constant_meas)
    b.free(slope);b.free(curve)
    b.mark('denormalize');normalize(b,x,flags,True);b.free(flags)
    return b.finish_circuit(),dict(events=b.events,peaks=dict(b.peaks))

def main():
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=('split','rows'),default='split')
    p.add_argument('--limit',type=int,default=84);p.add_argument('--keep',type=int,default=0);a=p.parse_args()
    tag=f'{a.mode}_q{a.limit}_k{a.keep}';path=HERE/'candidates'/f'{tag}.kmb'
    c,profile=build(**vars(a));c=compact(c,path);c.to_file(str(path),format='kmb')
    stats=c.operation_statistics();margin=json.loads((HERE/'results/model.json').read_text())['min_margin_units']
    nrot=stats.num_z_pow(min_clifford_hierarchy_level=4);tbase=stats.num_t_like
    synth=synthesis(nrot,margin,tbase)
    report=dict(status='BUILT',scope='complete circuit resource estimate; native verification pending',
        parameters=vars(a),Q=c.num_qubits,static_Toffolis=stats.num_toffoli_like,t_base=tbase,
        **synth,estimated_cost=c.num_qubits*(stats.num_toffoli_like-11+synth['synthesis_work']),
        profile=profile,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),path=str(path))
    (HERE/'results'/f'{tag}_build.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
