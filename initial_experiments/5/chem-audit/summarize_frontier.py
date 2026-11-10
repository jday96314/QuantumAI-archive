"""Checkpoint frontier table/plot; no claim of global nonexistence."""
import csv,json
from audit_core import HERE,write_json

rows=[];sources=[]
for path in sorted((HERE/'results').glob('scratch_*.json'))+sorted((HERE/'results').glob('robust_*.json')):
    data=json.loads(path.read_text())
    if data.get('rows'):sources.append((path.stem,data))
    for entry in data.get('rows',[]):
        for kind,metric in entry['metrics'].items():
            rows.append(dict(run=path.stem,stage=kind,step=entry['step'],**metric))
if rows:
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/'results/frontier_checkpoints.csv').open('w') as f:
        out=csv.DictWriter(f,fieldnames=keys);out.writeheader();out.writerows(rows)
    best=[]
    for name in sorted({r['run'] for r in rows}):
        rr=[r for r in rows if r['run']==name and r['stage']=='float']
        low=min(rr,key=lambda x:x['error_2'])
        bounded=[r for r in rr if r['Lambda']<=320]
        dominate=[r for r in rr if r['error_2']<=.02 and r['Lambda']<303.5283147431962]
        valid=[r for r in dominate if r['E_gap']>3.914898 and r['Lambda']>49 and r.get('gap_in_domain',True)]
        best.append(dict(run=name,best_float=low,best_bounded=min(bounded,key=lambda x:x['error_2']) if bounded else None,
            dominating_checkpoints=len(dominate),physically_dominating_checkpoints=len(valid)))
    scratch=[(name,data) for name,data in sources if name.startswith('scratch_')]
    guard=[r.get('training_norm_discrepancy') for _,data in sources for r in data['rows'] if r.get('training_norm_discrepancy') is not None]
    physical=[r.get('physical_contraction_discrepancy') for _,data in sources for r in data['rows'] if r.get('physical_contraction_discrepancy') is not None]
    write_json(HERE/'results/frontier_summary.json',dict(runs=best,number_of_checkpoints=len(rows)//3,metric_rows=len(rows),
        scratch_runs=len(scratch),expected_scratch_runs=18,all_scratch_complete=len(scratch)==18 and all(data['finished'] for _,data in scratch),
        scratch_dominating_checkpoints=sum(r['dominating_checkpoints'] for r in best if r['run'].startswith('scratch_')),
        maximum_recorded_training_norm_discrepancy=max(guard) if guard else None,
        maximum_recorded_physical_contraction_discrepancy=max(physical) if physical else None))
    import subprocess,sys
    subprocess.run([sys.executable,str(HERE/'plot_frontier.py')],check=True)
    print(json.dumps(best),flush=True)
