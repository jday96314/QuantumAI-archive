"""Measure compacted whole-circuit peaks and gates for workspace schedules."""
import argparse,itertools,json,time
from pathlib import Path
from build import build,save,HERE

p=argparse.ArgumentParser();p.add_argument('model',type=Path);p.add_argument('--wide',action='store_true');a=p.parse_args()
model=json.loads(a.model.read_text());cert=json.loads((HERE/'results'/f'{a.model.stem}_certificate.json').read_text())
model.update(affine_residual_range=cert.get('affine_residual_range'))
best=1e99;start=time.time();count=0;fails=0
if a.wide:
    grid=itertools.product([84,86,88,90,92,94],[0,2,4],[0,2,4],[14,17],[3,5])
else:
    grid=itertools.product([84,86,88,90,94],[0,2,4],[2],[17],[3,5])
with (HERE/'results'/f'{a.model.stem}_schedules.jsonl').open('w') as log:
 for limit,keep,delay,product,deferred in grid:
    if keep>=max(int(l[7]).bit_length() for l in model['leaves']):continue
    kw=dict(limit=limit,keep=keep,delay=delay,product_limit=product,deferred=deferred)
    try:c,profile=build(model,**kw)
    except (ValueError,AssertionError) as ex:
        fails+=1;log.write(json.dumps(dict(parameters=kw,status='INFEASIBLE',error=str(ex)))+'\n');log.flush();continue
    c,path,resources=save(c,a.model.stem+'_schedule_tmp',profile['stage_live_peaks'])
    t=c.operation_statistics().num_toffoli_like;q=c.num_qubits;cost=q*(t+cert['synthesis_toffoli_equivalent'])
    row=dict(Q=q,T_static=t,conservative_cost=cost,**resources,**profile)
    log.write(json.dumps(row)+'\n');log.flush();count+=1
    if cost<best:
        best=cost;target=HERE/'candidates'/f'{a.model.stem}_best.kmb';c.to_file(str(target),format='kmb')
        row.update(path=str(target),model=str(a.model),certificate=cert)
        (HERE/'results'/f'{a.model.stem}_best_build.json').write_text(json.dumps(row,indent=2)+'\n')
        print(json.dumps(dict(model=a.model.stem,Q=q,T_static=t,conservative_cost=cost,parameters=kw,seconds=time.time()-start)),flush=True)
print(json.dumps(dict(status='DONE',built=count,infeasible=fails,seconds=time.time()-start)),flush=True)
