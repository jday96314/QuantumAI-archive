"""Whole-oracle screens using the unchanged host confidence calculation.

These are NOT scores or native benchmark verification. All projected error
counts are explicitly labeled hypothetical, based on independent model runs.
"""
from pathlib import Path
import json,math
from dqi.verify_submission_dqi import clopper_pearson_upper_bound,repetition_overhead

HERE=Path(__file__).resolve().parent
BASE_Q=1004;BASE_T=1195802.1546875;BASE_COST=1435930717449.9397
CORE_T=824740;OLD_SCAN=371150;NATIVE_SHOTS=102400


def main():
    R0=repetition_overhead(clopper_pearson_upper_bound(0,NATIVE_SHOTS))
    assert math.isclose(R0,1.1912604648115832,rel_tol=1e-14)
    assert math.isclose(BASE_COST,R0*BASE_Q**2*BASE_T,rel_tol=1e-14)
    cf=json.loads((HERE/'results/cfraction_n70_profile.json').read_text())
    h=json.loads((HERE/'results/h_degree35_profile.json').read_text())
    small=json.loads((HERE/'results/cfraction_n6_m63_core-phase-cfraction-c.json').read_text())
    p=json.loads((HERE/'results/padding_pair_n35.json').read_text())
    pair=next(x['T'] for x in p if x['which']=='forward')
    cfQ=max(cf['Q'],840+small['complete_Q']-72)
    cfT=cf['T']+80*4095 # this emitted Forney scan does not skip zero masks
    # 4 pair adapters: A/B monomial -> squarefree, H/B -> monomial,
    # and both inverses. Reuse the incumbent's total 10,790T unpad budget
    # for H/B rotations and degree/control preparation. Assembly is pending.
    unpad=OLD_SCAN-88*4095
    hT=CORE_T+4*pair+h['T']+unpad+46*4095
    def screen(name,Q,T,errors=0):
        eps=clopper_pearson_upper_bound(errors,NATIVE_SHOTS);R=repetition_overhead(eps);cost=R*Q**2*T
        return dict(name=name,Q_assumption=Q,T_estimate=T,hypothetical_errors_in_102400=errors,
                    hypothetical_epsilon_plus=eps,hypothetical_R=R,cost_estimate=cost,
                    ratio_to_gap=cost/BASE_COST,conditional_speedup=BASE_COST/cost,
                    scope='resource/model screen, NOT a native result or task-cost improvement')
    screens=[screen('C-fraction, perfect chart handling at measured generic cost',cfQ,cfT),
             screen('Wronskian, perfect chart handling at measured generic cost',1004,hT),
             screen('Wronskian width sensitivity',1020,hT)]
    censuses=[]
    for label in ['a','b']:
        c=json.loads((HERE/f'results/census_complete_{label}.json').read_text());censuses.append(c)
        for kind,Q,T,key in [('C-fraction',cfQ,cfT,'cf_random_mask_phase_errors'),('Wronskian',1004,hT,'h_random_mask_phase_errors')]:
            errors=round(c[key]/c['shots']*NATIVE_SHOTS)
            r=screen(f'{kind}: hypothetical count extrapolated from model stream {label}',Q,T,errors)
            r.update(model_seed=c['seed'],model_shots=c['shots'],model_phase_errors=c[key],
                     caveat='Model extrapolation only; CP is evaluated at the hypothetical count, not claimed as observed. H assumes exact core cofactors.')
            screens.append(r)
    # An explicit extra H bank has a storage-policy bound, not a global bound.
    transport_oneway=sum(84*k+155 for k in range(1,35))
    retained_H_T=CORE_T+unpad+2*transport_oneway+46*4095
    screens.append(screen('carry a third H bank, favorable Q floor',840+420,retained_H_T))
    # Favorable direct-moment rolling circuit: share inverses across all
    # steps, count only 35 actual products/update, and omit every adapter.
    moment_lifetime=2*(4095-70)*35*42
    screens.append(screen('direct moment functional, favorable existing-70-moment bound',884,moment_lifetime))
    targets=[]
    for speedup in [1.1,1.2,1.5,2.0]:
        maxT=BASE_COST/(speedup*R0*cfQ**2)
        # Synthesis cost cannot be omitted; extraction would need this budget.
        synth=sum(s['T'] for s in cf['stages'] if s['name'].startswith('synthesize'))
        targets.append(dict(speedup=speedup,Q=cfQ,max_total_T=maxT,
                            max_oneway_CF_extraction_T=(maxT-80*4095-synth)/2,
                            current_oneway_CF_extraction_T=cf['stages'][0]['T']))
    sensitivities=[]
    for name,Q,T in [('C-fraction',cfQ,cfT),('Wronskian',1004,hT)]:
        k=0
        while repetition_overhead(clopper_pearson_upper_bound(k,NATIVE_SHOTS))*Q*Q*T<BASE_COST:k+=1
        sensitivities.append(dict(name=name,largest_hypothetical_error_count_that_beats_gap=k-1,
                                  first_noncompetitive_count=k))
    result=dict(status='NO_FULL_CANDIDATE_PROMOTED',baseline=dict(Q=BASE_Q,mean_T=BASE_T,R=R0,cost=BASE_COST),
                cf_core=cf,cf_scan_peak_extra_Q=44,cf_estimated_complete_Q=cfQ,cf_estimated_complete_T=cfT,
                H_charged_assembly=dict(core_T=CORE_T,pair_adapter_per_direction_T=pair,adapter_calls=4,
                    H_build_and_erase_T=h['T'],inherited_unpadding_budget_T=unpad,point_T=46,points=4095,
                    estimated_complete_T=hT,Q_assumption=1004,
                    limitation='H bank routing, gap dispatch, and composition with gap-core flags not integrated; no new cost assigned zero except unbuilt gap repair in the explicitly favorable screen'),
                screens=screens,CF_target_budgets=targets,confidence_sensitivity=sensitivities,
                direct_moment_field_products_only_T=moment_lifetime,
                H_recover_every_degree_oneway_T=sum(126*t*t+832*t-518 for t in range(1,18))+sum(126*t*t+664*t-582 for t in range(1,18)),
                acceptance='No complete full-size candidate passes or claims 102400-shot acceptance; incumbent unchanged.')
    (HERE/'results/resource_model.json').write_text(json.dumps(result,indent=2)+'\n')
    for s in screens:print(json.dumps(s),flush=True)
    print(json.dumps(dict(confidence_sensitivity=sensitivities,CF_target_budgets=targets)),flush=True)


if __name__=='__main__':main()
