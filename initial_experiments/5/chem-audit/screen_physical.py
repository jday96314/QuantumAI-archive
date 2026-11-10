"""List all saved integer states that warrant one-body/circuit preflight.

Reported pair-space/gap gates are only a shortlist. Native checks and the
actual chemistry cost decide selection; Lambda is not arbitrarily capped.
"""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
rows=[];seen=set()
def add(path,metrics,source):
    path=Path(path)
    if not path.exists() or str(path) in seen:return
    if metrics['error_2']>.02:return
    seen.add(str(path));rows.append(dict(fit=str(path.relative_to(HERE)),source=source,metrics=metrics,
        offset_repair_needed=metrics['E_gap']<=3.914898))
for path in sorted((HERE/'results').glob('scratch_*.json'))+sorted((HERE/'results').glob('robust_*.json')):
    data=json.loads(path.read_text());tag=path.stem
    if data.get('rows'):
        for kind,suffix in [('rounded','rounded'),('refit','integer')]:
            add(HERE/f'parameters/{tag}_{suffix}.npz',data['rows'][-1]['metrics'][kind],f'{tag}:last:{kind}')
            eligible=[r['metrics'][kind] for r in data['rows'] if r['metrics'][kind]['error_2']<=.02 and r['metrics'][kind]['E_gap']>3.914898 and r['metrics'][kind]['lambda_eff'] is not None]
            if eligible:
                best=min(eligible,key=lambda r:r['lambda_eff'])
                add(HERE/f'parameters/{tag}_best_physical_{kind}.npz',best,f'{tag}:best:{kind}')
        best_refit=min(data['rows'],key=lambda r:r['metrics']['refit']['error_2'])['metrics']['refit']
        add(HERE/f'parameters/{tag}_best_integer.npz',best_refit,f'{tag}:best_error_refit')
        passing_float=[r['metrics']['float'] for r in data['rows'] if r['metrics']['float']['error_2']<.02]
        if passing_float:
            best_float=min(passing_float,key=lambda r:r['Lambda'])
            add(HERE/f'parameters/{tag}_passing_float.npz',best_float,f'{tag}:passing_float_requires_round_check')
    elif 'metrics' in data and 'final' in data:
        add(HERE/f'parameters/{tag}.npz',data['metrics']['float'],f'{tag}:integer_refine')
    elif 'after' in data and 'arguments' in data:
        add(HERE/data['arguments']['out'],data['after']['rounded'],f'{tag}:gauged_rounding')
growth=HERE/'results/integer_growth.json'
if growth.exists():
    for r in json.loads(growth.read_text())['rows']:
        add(HERE/f'parameters/integer_growth_B{r["B"]}.npz',r['metrics']['float'],'integer_growth')
(HERE/'results/physical_shortlist.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
