import json
from pathlib import Path
from dqi.verify_submission_dqi import clopper_pearson_upper_bound,repetition_overhead
B=json.loads(Path('results/baseline_reproduction.json').read_text());rows=[]
for k in [0,1,2,4,8,380]:
 e=clopper_pearson_upper_bound(k,102400);r=repetition_overhead(e);cost=r*B['num_qubits']**2*B['toffoli_count'];rows.append(dict(phase_errors=k,eps_plus=e,R=r,cost_at_unchanged_Q_T=cost,gain_at_unchanged_Q_T=B['score']/cost,max_T_at_Q980_for_1p3_gain=B['score']/1.3/r/980**2))
Path('results/resource_budget.json').write_text(json.dumps(dict(scope='hypothetical budgets only; NOT replacement scores; unchanged baseline gates and width with substituted failure counts',rows=rows),indent=2)+'\n')
print(json.dumps(rows))
