"""Refresh the measured baseline table without touching other directories."""
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
lines=['## Measured scratch baseline checkpoints','',
    'Latest/terminal checkpoints at lambda_reg=300; errors are mHa. These are',
    'numerical two-body fits, not native passes. Other physical gates still apply.','',
    '| Shape | Step | Float error | Rounded error | Refit error | Float Lambda | Float gap | Refit Lambda |',
    '|---|---:|---:|---:|---:|---:|---:|---:|']
for R,B,C in [(1,719,350),(1,639,350),(2,447,175),(2,463,175),(2,479,175),(4,180,88)]:
    p=HERE/f'results/scratch_r{R}b{B}_l300.json'
    if not p.exists():continue
    z=json.loads(p.read_text());r=z['rows'][-1];f=r['metrics']['float'];q=r['metrics']['rounded'];s=r['metrics']['refit']
    lines.append(f'| R{R}/B{B}/C{C} | {r["step"]} | {1000*f["error_2"]:.4f} | {1000*q["error_2"]:.4f} | {1000*s["error_2"]:.4f} | {f["Lambda"]:.4f} | {f["E_gap"]:.6f} | {s["Lambda"]:.4f} |')
lines+=['','Reproduction: the first R1/B719 command is recorded below; all remaining',
    'baseline and regularization jobs are in `run_scratch.sh`, which skips',
    'completed baseline runs. `results/code/` preserves trainer/helper versions',
    'for the later runs. JSON arguments and NPZ hashes identify measured states.','',
    '```bash',
    'AUDIT_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 "$AUDIT_PY" frontier_fit.py --R 1 --B 719 --C 350 --seed 11 --lambda-reg 300 --weight-cap 5 --steps 20000 --checkpoint 2000 --precision float64 --checkpoint-device same --tag scratch_r1b719_l300',
    'bash run_scratch.sh',
    '# Prepared but unrun under the Phase B stop rule: run_robust.sh',
    '"$AUDIT_PY" summarize_frontier.py',
    '"$AUDIT_PY" screen_physical.py',
    '```','']
block='<!-- baseline-start -->\n'+'\n'.join(lines)+'\n<!-- baseline-end -->'
p=HERE/'RESULTS.md';text=p.read_text()
if '<!-- baseline-start -->' in text:
    start=text.index('<!-- baseline-start -->');end=text.index('<!-- baseline-end -->')+len('<!-- baseline-end -->')
    text=text[:start]+block+text[end:]
else:text+='\n\n'+block+'\n'
p.write_text(text)
