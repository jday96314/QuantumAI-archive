"""Gate-count screen for the emitted common matrix path, not a full score.

The proposed path is deliberately undercharged: only k-1 low coefficient
pairs are updated and coefficients are generated once per block. Upper bank
boundaries, completion, carries, and all serial repair are excluded. Compare
this partial cost to the incumbent's entire cofactor passes, including padding.
"""
import json
from pathlib import Path
from contextlib import contextmanager
from import_published import load
from low_workspace_inverse import bind
from blocked_monic_euclid import cofactor,schedule
from continued_fraction_series import apply as cf
from deferred_euclid_padding import fake_undo_ops,common_padding_ops
from controlled_linear import apply as linear
from mobius_series import BoundedArithmetic
from lookahead_six_quantum import apply as six
from matrix_transfer import apply as matrix

def setup():
    m=load();m.inv_word=bind(m);b=m.Builder()
    class Counter(m.Em):
        T=0
        def emit_mac(self,rt):self.T+=sum(term[0]=='P' for term in rt)
        def emit_free(self,rt):pass
    return m,b,Counter

def incumbent(k):
    m,b,Counter=setup();words=[b.qubits(12) for _ in range(2*k+6)]
    bp=words[:k-1];bb=words[k-1:2*k-1];q=words[2*k-1:2*k+5];pi=words[-1]
    flags=b.qubits(3);E=Counter(b,*b.qubits(3),m.Alloc(b));ff,_=schedule(bp,bb,q)
    passes=[]
    for back,norm in [(False,False),(True,False),(False,True)]:
        before=E.T+b._cur_stats['toffoli']
        for j in (reversed(range(3)) if back else range(3)):
            pp,bb,be,al=ff[j];ops=fake_undo_ops(pp,bb,be,al)+common_padding_ops(pp,bb,[be,al],False)
            if back:linear(E.b,ops,flags[j],True)
            cofactor(m,E,pp,bb,be,al,pi,back,norm)
            if not back:linear(E.b,ops,flags[j])
        passes.append(E.T+b._cur_stats['toffoli']-before)
    return sum(passes),passes

def proposed(n):
    m,b,Counter=setup();rows=[[b.qubits(12) for _ in range(k)] for k in (2,3,3,3)]
    pp=[b.qubits(12) for _ in range(n)];bb=[b.qubits(12) for _ in range(n)]
    E=Counter(b,*b.qubits(3),m.Alloc(b));matrix(m,E,rows,pp,bb)
    return E.T+b._cur_stats['toffoli'],b.next_qubit

def head(which,record=True):
    m,b,Counter=setup();words=[b.qubits(12) for _ in range(6)];flags=b.qubits(3)
    E=Counter(b,*b.qubits(3),m.Alloc(b));math=BoundedArithmetic(m,E);math.virtual_zero=True
    if which=='cf':cf(m,E,math,words,flags=flags,record=record)
    else:six(m,E,*words)
    return E.T+b._cur_stats['toffoli']

def gap_branch(k):
    from gap_block import apply as gap
    m,b,Counter=setup();words=[b.qubits(12) for _ in range(70)]
    E=Counter(b,*b.qubits(3),m.Alloc(b));math=BoundedArithmetic(m,E);math.virtual_zero=True
    gap(m,E,math,words[:k-1],words[k-1:2*k-1],words[2*k-1:-1],words[-1])
    return E.T+b._cur_stats['toffoli']

def main():
    rows=[]
    for k in range(1,22,3):
        base,passes=incumbent(k);new,q=proposed(k-1)
        rows.append(dict(k=k,interior_words_per_bank=k-1,baseline_all_cofactor_passes=base,baseline_passes=passes,new_one_interior_pass=new,new_three_interior_passes=3*new,delta_forward=3*new-base,component_Q=q))
    delta=2*sum(r['delta_forward'] for r in rows)
    cfT,cfNoRecord,sixT=head('cf'),head('cf',False),head('six')
    head_delta=2*len(rows)*(2*sixT-cfT-cfNoRecord)
    bridge=json.loads(Path('results/matrix_coefficients.json').read_text())
    assert bridge['status']=='PASS'
    bridgeT=bridge['circuits'][0]['T'];generation=2*len(rows)*bridgeT
    partial=delta+head_delta+generation
    report=dict(scope=__doc__,status='KILL_THIS_CONSTRUCTION' if partial>90000 else 'CONTINUE',rows=rows,forward_inverse_interior_delta=delta,baseline_head_cf=cfT,baseline_head_cf_without_record=cfNoRecord,adjacent_minor_head=sixT,head_delta_estimate=head_delta,coefficient_generation_measured_T=bridgeT,one_generation_per_block_delta=generation,partial_projected_delta=partial,partial_projected_mean_T=1019865.856484375+partial,including_one_coefficient_cleanup_per_block_delta=partial+generation,including_one_coefficient_cleanup_projected_mean_T=1019865.856484375+partial+generation,bridge_Q=bridge['circuits'][0]['Q'],optimistic_embedded_bridge_Q=840+bridge['circuits'][0]['Q']-72+6+3,excluded=['top boundary coefficient updates','additional coefficient work required across cofactor passes','completion handling','carry alignment and cleanup','serial exception repair','degree-three and (1,2) local head support'],not_a_lower_bound_on_other_architectures=True)
    Path('results/transfer_resource_screen.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
    branches=[dict(k=k,uncontrolled_gap_block_T=gap_branch(k)) for k in range(1,22,3)]
    selector=dict(scope='Screen of naive quantum selection between two complete arithmetic blocks. All arithmetic in both branches is emitted. Costs for controlling gates, completion, final slots, and serial repair are omitted. This is not a lower bound on a selector that shares arithmetic.',rows=branches,optimistic_added_T=2*sum(r['uncontrolled_gap_block_T'] for r in branches),status='KILL_NAIVE_WHOLE_BLOCK_SELECTOR')
    Path('results/whole_block_selector_screen.json').write_text(json.dumps(selector,indent=2)+'\n')

if __name__=='__main__':main()
