"""Unmodified host CP objective, conditional census and exact splice screens."""
import json, math
from pathlib import Path
from dqi.verify_submission_dqi import clopper_pearson_upper_bound, repetition_overhead
ROOT=Path(__file__).resolve().parent
C0=1598204229384.9265
R0=repetition_overhead(clopper_pearson_upper_bound(0,102400))
RR=repetition_overhead(clopper_pearson_upper_bound(380,102400))
P=4025*4095/4096

def row(name,q,t,scope):
    c=R0*q*q*t
    return dict(name=name,Q=q,T=t,R0=R0,cost_R0=c,ratio_R0=C0/c,Rref=RR,cost_Rref=RR*q*q*t,ratio_Rref=C0/(RR*q*q*t),scope=scope)

def main():
    profiles=json.loads((ROOT/'results/profiles.json').read_text())
    F=next(r for r in profiles if r['kind']=='fprime' and r['d']==35)
    sq=next(r for r in profiles if r['kind']=='square' and r['d']==35)
    # Exact decomposition of the SAME deterministic emitted Frobenius algorithm.
    mulT=23352
    ix=2*(11*sq['T']+10*mulT)
    assert F['T']==2*(ix+2*mulT)
    # Width: public A,B,retained work=3d; derivative,c=2d;
    # all21 inversion intermediates; convolution=2d-1; emitter pivots=3.
    assert F['Q']==12*(3*35+2*35+21*35+69)+3
    # Public data is already 105words vs incumbent70; use component peak as
    # a favorable lower bound on complete decoder width, omit flags entirely.
    tested=[row('A_Frobenius_complete_lifetime_favorable_embedding',F['Q'],649802+9792+11*P+2*F['T'],'exact component charges; decoder-control embedding and new gap repair omitted'),
            row('C_70_full_width_reference_updates',1707,2*70*31962+80*P+10000,'deterministic uniform disjoint update schedule; additional correction computation/branch history omitted; does not claim recurrence correctness'),
            row('A_no_conversion_lower_screen',961,649802+9792+11*P,'conversion=0 and gap repair=0 used only for favorable headroom'),
            row('wide_BM_Psi',1837,442397+23*(P-4024),'published fixed-mask T with expected point-count correction; not a host-verified circuit'),
            row('H_generic_Wronskian',961,649802+2*(7140+62098)+46*P,'generic adapter; H division, total singular/gap handling omitted'),
            row('D_general_displacement',2200,2*42*3*35**2+2*35*66+20000+80*P,'structured-elimination arithmetic estimate; all four generator vectors counted'),
            row('F_indexed_root_list',2200,300000+35*4094,'optimistic decoding/factoring allowance, lookup cleanup may add gates'),
            row('G_streaming_residue',1360,649802+2*42*4095*35,'no quotient history or reduction adapters; favorable screen')]
    targets=[row('C_resident_rational_planning',961,280000+10000+80*P,'unimplemented core allowance; reference update does not attain it'),
             row('C_2p5_target',940,240000+10000+80*P,'unimplemented core/width target; detailed build plan is conditional'),
             row('A_direct_target',961,649802+9792+11*P+200000,'100k per converter direction; G=0, not a measured estimate'),
             row('A_trimmed_counter_target',948,624206+9792+11*P+60000,'30k converter per direction, assumed flags and G=0; unimplemented'),
             row('D_compressed_generator_target',1100,2*42*3*35**2+2*35*66+20000+80*P,'unimplemented compression; no claim symmetry removes diagonal data')]
    handoff=json.loads((ROOT/'results/lazy_handoff_screen.json').read_text())
    for r in handoff:
        r.update(row(f"B_handoff_k{r['k']}",r['Q_lower_bound'],r['total_T_zero_cost_handoff'],r['scope']))
    targets.sort(key=lambda r:r['cost_R0']);tested.sort(key=lambda r:r['cost_R0']);handoff.sort(key=lambda r:r['cost_R0'])
    budgets=[]
    for q in [880,940,948,960,961,1050,1200,1383,1837]:
        total25=C0/(2.5*R0*q*q);total2=C0/(2*R0*q*q)
        budgets.append(dict(Q=q,T_limit_2x=total2,T_limit_2p5x=total25,C_core_plus_adapter_limit_2p5x=total25-80*P,A_converter_one_way_limit_2x=(total2-649802-9792-11*P)/2))
    sensitivity=[]
    for k in [0,1,8,15,50,380]:
        eps=clopper_pearson_upper_bound(k,102400);R=repetition_overhead(eps)
        sensitivity.append(dict(failures=k,eps_plus=eps,R=R,C_target940_cost=R*940**2*(240000+10000+80*P)))
    out=dict(status='NO_ACCEPTED_REPLACEMENT',R0=R0,Rref=RR,active_points=P,incumbent_cost=C0,initial_targets_ranked=targets,priced_realizations_ranked=tested,lazy_handoff_ranked=handoff,budgets=budgets,failure_sensitivity=sensitivity,Frobenius_exact_decomposition=dict(square_T=sq['T'],multiply_T=mulT,inverse_xor_T=ix,fprime_one_direction_T=F['T'],fprime_lifetime_T=2*F['T'],component_peak=F['Q']),qualification='Conditional objective screens, not full native DQI scores or posterior expectations. Native components do not establish decoder failure rates. No family dropped solely on a losing first draft.')
    (ROOT/'results/budget.json').write_text(json.dumps(out,indent=2)+'\n')
    for r in targets+handoff[:2]+tested:print(r['name'],r['Q'],round(r['T'],3),f"{r['cost_R0']:.6e}",f"{r['ratio_R0']:.4f}x")

if __name__=='__main__':main()
