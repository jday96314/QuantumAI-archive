"""Cross-check native components against the count screen; price projections."""
from collections import Counter
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
def read(name):return json.loads((HERE/'results'/name).read_text())

def main():
    base=read('attribution_b3_s22.json');opt=read('attribution_b3_s22_trimmed_retain.json')
    screen=read('schedule_screen.json');native=read('native_blocks.json');middle=read('native_middle.json')
    for row in native['rows']:
        kind='block' if row['block']==3 else 'serial'
        exact=next(r for r in screen['rows'] if r.get('k')==row['k'] and r['size']==row['block'] and r['kind']==kind)
        expected=exact['forward_T'] if row['mode']=='forward' else exact['reverse_T'] if row['mode']=='inverse' else exact['forward_T']+exact['reverse_T']
        assert row['T']==expected,(row,exact)
    core=base['forward_T']+base['reverse_T'];new=opt['forward_T']+opt['reverse_T'];scan=370152
    assert screen['optimal_core_T']==new
    arithmetic=Counter()
    for row in base['arithmetic_calls']:arithmetic[(row['region'].split('/')[0],row['area'],row['operation'])]+=row['T']
    projected=1.2*961**2*(new+scan)
    hypothetical=1.2*852**2*(new+scan)
    report=dict(status='COMPONENT_PASS_NOT_ACCEPTED_DECODER',baseline_core_T=core,prototype_core_T=new,
        core_saving_T=core-new,core_saving_fraction=(core-new)/core,selected_and_prototype_peak_Q=961,
        scan_static_T=scan,baseline_static_total_T=core+scan,prototype_static_total_T=new+scan,
        complete_static_saving_fraction=(core-new)/(core+scan),
        objective_projection_R1_2_Q961=projected,objective_projection_R1_2_Q852_hypothetical=hypothetical,
        objective_projection_assumption='static total used as mean-T approximation; R=1.2 assumed, not measured; Q852 not achieved',
        native_component_circuits=len(native['rows'])+len(middle['rows']),shots_per_component=native['shots'],
        native_component_simulated_shots=native['shots']*len(native['rows'])+middle['shots']*len(middle['rows']),
        local_schedule_components=120,local_schedule_directions=240,
        gu_61_percent_reverse_hypothetical_core_T=324901*1.61,
        gu_61_percent_reverse_hypothetical_saving_fraction=.39/2,
        current_asymmetric_reverse_saving_estimate_T=0,
        arithmetic_operations=[dict(direction=k[0],area=k[1],operation=k[2],T=v) for k,v in arithmetic.items()],
        stopped='less than 10% reducible within tested compact-state arithmetic/schedule family; no global decoder lower bound',
        full_native_102400_and_independent_salts_run=False,selected_replaced=False)
    (HERE/'results/checkpoint.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
