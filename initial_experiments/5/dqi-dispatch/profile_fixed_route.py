"""Cost screen of the width-feasible separate carry-chart realization.

Counts the native shared slot and Schur head, regular controlled moments,
one carried eight-moment chart, and one carried ten-moment chart (nine after
removing a known leading zero). Each is emitted on every shot. The chart
normalizations, upper grafts, mode computation, serial corrections, completion
and gap-history cleanup remain excluded. This is not a bound on a common
Möbius path that shares arithmetic across carry charts.
"""
import json,math
from pathlib import Path
from import_published import load
from low_workspace_inverse import bind
from core_components import TrimmedArithmetic
from mobius_series import apply as ordinary_mobius
from shared_slot import apply as shared
from moment_probe import emit as moment

def setup():
    m=load();m.inv_word=bind(m);b=m.Builder()
    class Counter(m.Em):
        T=0
        def emit_mac(self,rt):self.T+=sum(t[0]=='P' for t in rt)
        def emit_free(self,rt):pass
    return m,b,Counter

def scale_slot(n,back,norm,routed):
    m,b,Counter=setup();X=[b.qubits(12) for _ in range(n)];Y=[b.qubits(12) for _ in range(n)];c=b.qubits(12);p=b.qubits(12);flag=b.qubits(1)[0]
    E=Counter(b,*b.qubits(3),m.Alloc(b))
    if routed:shared(m,E,X,Y,c,p,flag,back,norm)
    else:
        with m.inv_word(E,c) as ci,m.opened(E,'reference_scale_product',12) as t:
            for x in X:m.mon_slot(E,x,ci if back else c,c if back else ci,t)
            if norm:m.mon_slot(E,p,c if back else ci,ci if back else c,t)
    return E.T+b._cur_stats['toffoli']

def ordinary_head(k):
    m,b,Counter=setup();A=[b.qubits(12) for _ in range(k+1)];D=[b.qubits(12) for _ in range(k)];p=b.qubits(12);T=[b.qubits(12) for _ in range(6)]
    E=Counter(b,*b.qubits(3),m.Alloc(b));a=TrimmedArithmetic(m,E);a.virtual_zero=True
    ordinary_mobius(m,E,a,A,D,p,T,p_is_inverse=True)
    return E.T+b._cur_stats['toffoli']

def rotations(r):
    # Canonical partial prefix -> working A/B prefix. One direction only;
    # returning to the canonical frame is omitted to undercharge the screen.
    cases=[]
    # Modes 1,4,5 all move a trailing zero into the same missing target
    # coordinate. Modes 4/5 have two zeros; extending their rotation by the
    # second zero gives exactly the same valid chart. Emit this map once.
    if r>=2:cases.append(('1|4|5',r+1,1))
    if r>=3:cases.append((2,r+2,2))
    if r>=2:cases.append((6,r,2))
    m,b,Counter=setup();words=[b.qubits(12) for _ in range(2*r+2)];flag=b.qubits(1)[0];rows=[]
    for mode,n,amount in cases:
        before=b._cur_stats['toffoli'];seen=set()
        for start in range(n):
            if start in seen:continue
            cycle=[];i=start
            while i not in seen:seen.add(i);cycle.append(i);i=(i+amount)%n
            for i in reversed(cycle[1:]):
                for x,y in zip(words[start],words[i]):b.cswap(flag,x,y)
        T=b._cur_stats['toffoli']-before
        assert T==12*(n-math.gcd(n,amount))
        rows.append(dict(mode=mode,words=n,amount=amount,T=T))
    return rows

def main():
    heads=[];slots=[];routing=[]
    for k in range(1,22,3):
        base=ordinary_head(k);normal=moment(k,6,'forward',True,True,False)
        carry1=moment(k,8,'forward',True,True);carry2=moment(k,9,'forward',True,True)
        heads.append(dict(k=k,baseline_moment=base,conditional_regular=normal,carried_one=carry1,carried_two=carry2,
                          full_decoder_delta=4*(normal['T']-base+carry1['T']+carry2['T'])))
        for j in range(3):
            r=k+j;n=r-1;delta=0;details=[]
            for back,norm in [(False,False),(True,False),(False,True),(True,True),(False,False),(True,False)]:
                old,new=scale_slot(n,back,norm,False),scale_slot(n,back,norm,True);delta+=new-old;details.append(dict(inverse=back,norm=norm,old=old,new=new))
            slots.append(dict(k=r,n=n,delta=delta,passes=details))
            rs=rotations(r);routing.append(dict(k=r,one_way_per_pass=sum(x['T'] for x in rs),cases=rs))
    head_delta=sum(r['full_decoder_delta'] for r in heads)
    slot_delta=sum(r['delta'] for r in slots)
    # Six passes: three cofactor passes in each direction.
    route_delta=6*sum(r['one_way_per_pass'] for r in routing)
    schur=json.loads(Path('results/local_schur.json').read_text());assert schur['status']=='PASS'
    schurT=schur['circuits'][0]['T'];schur_delta=14*(2*schurT-969-936)
    added=-25596+head_delta+slot_delta+schur_delta+route_delta
    report=dict(status='KILL_SEPARATE_CHART_HEADS' if added>300000 else 'CONTINUE',scope=__doc__,heads=heads,slots=slots,routing=routing,
                reused_core_saving=25596,conditional_moment_delta=head_delta,shared_slot_delta=slot_delta,routed_schur_delta=schur_delta,
                one_way_bank_routing_delta=route_delta,partial_projected_added_T=added,partial_projected_mean_T=1019865.856484375+added,
                native_component_embedding_peak=max(x[n]['Q'] for x in heads for n in ['conditional_regular','carried_one','carried_two']),
                planned_full_peak=973,excludes=['carry normalization','chart return rotations','upper cofactor grafts','mode derivation','gap-history updates/cleanup','completion detection/padding','serial correction','pending head Schur changes'],
                accepted_full_decoder=False,not_a_lower_bound_on_a_common_moment_map=True)
    adopted=Path('parameters/base_adoption.json')
    if adopted.exists():
        base=json.loads(adopted.read_text());assert base['status']=='ADOPTED_READ_ONLY_ACCEPTED_BASE'
        report.update(base_accepted=True,base_circuit=base['circuit'],base_sha256=base['sha256'],
                      partial_added_T_over_accepted_base=added+25596,
                      partial_projected_mean_T=base['default_mean_T']+added+25596)
    Path('results/fixed_route_resource_screen.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['heads','slots','routing']}),flush=True)

if __name__=='__main__':main()
