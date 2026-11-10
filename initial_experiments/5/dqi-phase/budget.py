"""Objective screens using host Clopper-Pearson and measured components.

An optimistic screen is NOT an assembled circuit or native cost claim.
In particular, unknown compact-prior -> omega conversion is never called free
in a proposed accepted implementation: zero is used only to reject the route
even under a favorable lower-cost assumption.
"""
import json
import math
from pathlib import Path
import numpy as np
from dqi.verify_submission_dqi import clopper_pearson_upper_bound,repetition_overhead
from oracle import (EXP,mul,vmul,power,trace,interpolation_inverse,
                    MASK_TO_DUAL,DUAL_TO_MASK)
from tower_inverse import map_ops


def main():
    control=json.loads(Path('results/baseline_audit.json').read_text())
    psi=json.loads(Path('results/psi_schedule.json').read_text())
    fast=json.loads(Path('results/fast_modular.json').read_text())
    f35=next(r for r in fast['profiles'] if r['degree']==35)
    Q0=control['native_Q'];T0=control['native_mean_T'];C0=control['native_cost']
    eps=clopper_pearson_upper_bound(control['native_phase_errors']+control['native_value_errors'],102400)
    R=repetition_overhead(eps)
    assert math.isclose(eps,control['native_epsilon_plus'],rel_tol=1e-14)
    assert math.isclose(R,control['native_R'],rel_tol=1e-14)
    fixed=control['core_forward_and_reverse_T']+control['cofactor_unpadding_and_repadding_T']
    active=(4095-70)*(4095/4096)
    def screen(label,Q,T,scope):
        cost=R*Q*Q*T
        return dict(label=label,Q=Q,mean_T=T,R_assumption=R,cost_estimate=cost,cost_ratio_to_incumbent=cost/C0,
                    speedup_estimate=C0/cost,scope=scope,
                    absolute_R1_cost=Q*Q*T,
                    zero_failures_102400_CP_cost=repetition_overhead(clopper_pearson_upper_bound(0,102400))*Q*Q*T,
                    max_Q_for_break_even=math.sqrt(C0/(R*T)),max_Q_for_1_25x=math.sqrt(C0/(1.25*R*T)))
    screens=[screen('80T division, evaluator already present',Q0,fixed+80*4095*(4095/4096),'hypothetical favorable screen; wrong layout for current compact scan, adapter excluded'),
             screen('80T division plus per-point compact-prior inversion',Q0,fixed+182*4095*(4095/4096),'Toffoli screen; 102T inverse lifetime + 80T body, width optimistically unchanged'),
             screen('23T Psi, evaluator conversion assumed zero cost',psi['live_peak'],fixed+psi['roundtrip_T']+23*active,'optimistic rejection screen, not complete circuit; actual Psi liveness + full inverse; missing evaluator conversion can only add work'),
             screen('11T F, fast in-place omega transform',f35['F_lifetime_Q'],fixed+psi['roundtrip_T']+f35['full_F_overwrite_lifetime_T']+11*active,'optimistic rejection screen; Q is F component alone; excludes compact-prior conversion, variable-degree adapter and Psi degree-35 reduction'),
             screen('11T F, fast out-of-place omega transform',f35['Q'],fixed+psi['roundtrip_T']+f35['out_of_place_F_lifetime_T']+11*active,'optimistic rejection screen; Q is modular MAC component alone; retains original omega and new F bank; missing adapters excluded')]
    ideal=[]
    for body in [88,80,23,11]:
        T=fixed+body*active
        ideal.append(screen(f'{body}T point, no preparation at all',Q0,T,'unachievable best-case headroom bound, preparation deliberately zero'))
    # Count every CNOT in the implemented classical mask preparation/unwind
    # and conversion for the full 4095-position scan, without serializing it.
    inv=interpolation_inverse();pc=np.array([sum(mul(c,1<<i).bit_count() for i in range(12)) for c in range(4096)])
    powers=np.array([[power(2,(i+1)*j) for j in range(4095)] for i in range(70)])
    dual=len(map_ops([int(MASK_TO_DUAL[1<<i]) for i in range(12)]))
    coord=len(map_ops([int(DUAL_TO_MASK[1<<i]) for i in range(12)]))
    classical=dict(extra_quantum_qubits=0,extra_classical_bits=840,Toffolis=0,
        guaranteed_skips=70,expected_total_skips=70+4025/4096,
        interpolation_and_residual_CX=2*int(pc[inv].sum()+pc[powers].sum()),
        dual_basis_CX_or_SWAP=2*4095*dual,residual_coordinate_CX_or_SWAP=2*4025*coord,
        syndrome_Z_instructions=70*sum(trace(mul(1<<i,1<<j)) for i in range(12) for j in range(12)),
        note='CX/SWAP counts are instructions under the exact free linear maps; classical operations and Z gates have zero objective T/Q cost, but are emitted and reversed')
    # Sensitivity uses independent hypothetical failure counts; no hash search.
    sens=[]
    for k in [300,380,460,600]:
        e=clopper_pearson_upper_bound(k,102400);r=repetition_overhead(e)
        sens.append(dict(failures=k,shots=102400,epsilon_plus=e,R=r,control_QT_cost=r*Q0*Q0*T0))
    joint=json.loads(Path('../../4/codex/results/joint_cubic_trace_product_n70_m4095_verify.json').read_text())
    headroom=[]
    for Q in [961,1000,1050,1100]:
        limit=C0/(1.25*R*Q*Q)
        headroom.append(dict(Q=Q,max_total_mean_T_for_1_25x=limit,
            all_preparation_and_inverse_T_budget_at_23_point=limit-fixed-23*active,
            all_preparation_and_inverse_T_budget_at_11_point=limit-fixed-11*active,
            remaining_F_and_adapters_after_measured_Psi=limit-fixed-11*active-psi['roundtrip_T']))
    out=dict(status='NO_FULL_CANDIDATE_PROMOTED',control=control,host_CP_reproduction=dict(epsilon_plus=eps,R=R),
             fixed_non_point_T=fixed,expected_retained_active_positions=active,
             component_screens=screens,ideal_headroom=ideal,preparation_headroom=headroom,classical_mask_overhead=classical,failure_sensitivity=sens,
             joint_87T_prior_native=dict(shots=joint['verified_shots'],Q=joint['num_qubits'],mean_T=joint['toffoli_count'],R=joint['repetition_factor'],cost=joint['score'],cost_ratio=joint['score']/C0),
             acceptance='No >=15-20% realistic improvement in screened schedules; no full integration or new 102400-shot acceptance claimed. Existing selected artifacts preserved.')
    Path('results/budget.json').write_text(json.dumps(out,indent=2)+'\n')
    for row in screens:print(json.dumps(row),flush=True)


if __name__=='__main__':main()
