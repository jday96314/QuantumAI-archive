"""Write the final measured findings and explicit conditional stopping scope."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
summary=json.loads((HERE/'results/frontier_summary.json').read_text())
assert summary['all_scratch_complete'] and summary['scratch_dominating_checkpoints']==0
selected=json.loads((HERE/'selection.json').read_text())
data=[]
for p in (HERE/'results').glob('scratch_*.json'):
    z=json.loads(p.read_text())
    if z.get('rows'):
        for r in z['rows']:data.append((p.stem,r))
lines=['## Final Phase B findings and stopping decision','',
    '**Stopped under the user-defined condition.** The 46.7861796 bound is',
    'universal for the native fixed-constant problem, and none of the tested',
    'from-scratch Lambda-regularized float fits reaches error2 <=20 mHa with',
    'Lambda below the incumbent 303.5283147. There are six random scratch',
    'starts (20,000 steps each) and twelve regularization continuations',
    '(6,000 steps each). The continuations share those six basins and follow',
    '20-sweep column-gauge reductions; they are not independent random starts.','',
    '| Shape | Best float error (mHa) | Lambda there | Best float error at Lambda <=320 | Lambda there | Rounded error there | Refit error there | Refit Lambda |',
    '|---|---:|---:|---:|---:|---:|---:|---:|']
for R,B,C in [(1,719,350),(1,639,350),(2,447,175),(2,463,175),(2,479,175),(4,180,88)]:
    rs=[r for name,r in data if name.startswith(f'scratch_r{R}b{B}_l')]
    best=min(rs,key=lambda r:r['metrics']['float']['error_2'])['metrics']['float']
    bounded=[r for r in rs if r['metrics']['float']['Lambda']<=320]
    low=min(bounded,key=lambda r:r['metrics']['float']['error_2'])['metrics']
    f,q,s=low['float'],low['rounded'],low['refit']
    lines.append(f'| R{R}/B{B}/C{C} | {best["error_2"]*1000:.4f} | {best["Lambda"]:.4f} | {f["error_2"]*1000:.4f} | {f["Lambda"]:.4f} | {q["error_2"]*1000:.4f} | {s["error_2"]*1000:.4f} | {s["Lambda"]:.4f} |')
minimum=min(r['metrics']['refit']['error_2'] for _,r in data)
lines+=['',
    f'The scratch-family minimum rounded/refit error over all saved checkpoints is **{minimum*1000:.4f} mHa**.',
    f'The complete CSV has **{summary["number_of_checkpoints"]} checkpoint triples / {summary["metric_rows"]} metric rows**.',
    'Every triple records float, direct rounding and exact rank-C refit error,',
    'Lambda, gap, effective normalization and weight norms. Saved research',
    'NPZs, JSON arguments and source snapshots identify the exact measured states.','',
    f'The maximum recorded float32-vs-float64 error-norm discrepancy is **{summary["maximum_recorded_training_norm_discrepancy"]*1000:.6f} mHa**;',
    f'the maximum recorded physical-contraction discrepancy is **{summary["maximum_recorded_physical_contraction_discrepancy"]:.3g} Ha**.',
    'These are numerical consistency checks, not native acceptance. The final',
    '`physical_shortlist.json` is empty: no scratch integer state meets the',
    'two-body accuracy gate, even before building one-body words and oracles.','',
    'The expected-half-LSB-noise, integer-descent and integer-growth branch was',
    '**prepared but not executed**, because the explicit Phase B stop condition',
    'was met after the frontier audit. `run_robust.sh`, `integer_refine.py` and',
    '`integer_grow.py` are unrun research implementations, not evidence of',
    'success or failure for those methods. The CPU preflight is an implementation',
    'check only and is not counted as a completed frontier fit.','',
    'The negative fit results apply to these seeds, budgets, parameterization',
    'and regularization paths. They do not establish global infeasibility of',
    'these shapes or a universal chemistry-cost bound. Indeed the preserved',
    'known integer basin is feasible at R1/B719; the random searches did not',
    'reproduce it. Other basins, gauges, longer optimization and oracle designs',
    'remain possible. The rigorous universal result concerns normalization only.','',
    f'**Final selection remains `candidates/banked`, cost {selected["score"]:,.4f}, improvement {selected["improvement"]:.9f}x.**',
    'It passes two full-default native salts with 50,880 verified shots each.',
    'The <=567M target and <=420M stretch were not reached. No later research',
    'fit replaces the native-verified margin gain.','',
    'Artifacts: [all checkpoint values](results/frontier_checkpoints.csv),',
    '[frontier plot](results/frontier.png), [vector plot](results/frontier.svg),',
    '[native result](results/banked_native.json), [independent native result](results/banked_native_independent.json),',
    '[selection](selection.json), [source/artifact hashes](artifact_manifest.json).','']
block='<!-- final-findings-start -->\n'+'\n'.join(lines)+'\n<!-- final-findings-end -->'
p=HERE/'RESULTS.md';text=p.read_text()
if '<!-- final-findings-start -->' in text:
    start=text.index('<!-- final-findings-start -->');end=text.index('<!-- final-findings-end -->')+len('<!-- final-findings-end -->')
    text=text[:start]+block+text[end:]
else:text+='\n\n'+block+'\n'
text=text.replace('## Phase B implementation and ongoing checkpoints','## Phase B implementation and finite search')
status='**Phase B stopped after the requested frontier audit.** The bound is universal; no tested scratch fit beats the incumbent error/Lambda point. The prepared noise/descent/growth branch remains unrun under this stopping rule.\n\n'
if status not in text:text=text.replace('## Hypotheses and protocol',status+'## Hypotheses and protocol')
p.write_text(text)
