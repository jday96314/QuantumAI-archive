"""Exact resource search over variable one/two/three-slot schedules."""
import json,multiprocessing as mp
from pathlib import Path
from import_published import load
from low_workspace_inverse import bind
from mobius_series import BoundedArithmetic
from profile_hybrid import DiscardTape
from gap_serial import choose_wait
from fused_gap_block import apply as fused,previous_shift
from onepass_gap import apply as onepass
from dqi.verify_submission_dqi import clopper_pearson_upper_bound,repetition_overhead


def stage(task):
    k,count=task;m=load();m.inv_word=bind(m);b=m.Builder();b._items=DiscardTape();W=[b.qubits(12) for _ in range(70)];Z=[(),b.qubits(12)];terminal,phase=b.qubits(2)
    waits=b.qubits(k+count);flags=b.qubits(k+count);E=m.Em(b,*b.qubits(3),m.Alloc(b));math=BoundedArithmetic(m,E);math.virtual_zero=True
    bp=W[:k-1];bb=W[k-1:2*k-1];c=W[2*k-1:-1];p=W[-1]
    if count==1:
        with previous_shift(m,E,k,k,waits,flags,terminal) as g:
            choose_wait(m,E,c[-1],c[-2],waits[k-1],g,waits[k])
            onepass(m,E,math,bp,bb,c,p,Z,waits[k-1],waits[k],g,flags[k])
    else:fused(m,E,math,bp,bb,c,p,Z,k,count,waits,flags,terminal)
    for i in range(k,k+count):
        if i>=2:E.b.ccx(waits[i-2],flags[i],terminal)
    return dict(k=k,count=count,Q=b.peak_qubits,T_forward=b._cur_stats['toffoli'],T=2*b._cur_stats['toffoli'])


if __name__=='__main__':
    cache=Path('results/schedule_stage_profiles.json');rows=json.loads(cache.read_text()) if cache.exists() else []
    done={(r['k'],r['count']) for r in rows};tasks=[(k,c) for k in range(1,35) for c in [1,2,3] if k+c<=35 and (k,c) not in done]
    with mp.Pool(3) as pool:
        for row in pool.imap_unordered(stage,tasks):
            rows.append(row);cache.write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(row),flush=True)
    costs={(r['k'],r['count']):r for r in rows}
    ref=json.loads(Path('results/profile_b1_after1_laurent.json').read_text());fixed=ref['T']-sum(r['T'] for r in ref['rows'] if r['label'].startswith(('forward_','inverse_')))
    scan=json.loads(Path('results/scan_cost.json').read_text())['scan_T']+182
    R=repetition_overhead(clopper_pearson_upper_bound(0,102400));out=[]
    for cap in sorted({975,*[r['Q'] for r in rows if r['Q']>=975]}):
        dp={35:(0,[])}
        for k in reversed(range(1,35)):
            options=[]
            for count in [1,2,3]:
                r=costs.get((k,count))
                if r and r['Q']<=cap and k+count in dp:options.append((r['T']+dp[k+count][0],[count]+dp[k+count][1]))
            if options:dp[k]=min(options)
        if 1 not in dp:continue
        T=dp[1][0]+fixed+scan;cost=R*cap*cap*T
        out.append(dict(Q=cap,T=T,core_T=dp[1][0]+fixed,scan_T=scan,R=R,cost=cost,gain=1598204229384.9265/cost,counts=dp[1][1]))
    out.sort(key=lambda r:r['cost']);Path('results/schedule_screen.json').write_text(json.dumps(out,indent=2)+'\n')
    Path('parameters/best_schedule.json').write_text(json.dumps(out[0],indent=2)+'\n');print(json.dumps(out[0]),flush=True)
