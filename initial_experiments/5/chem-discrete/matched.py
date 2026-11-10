"""Serial matched exact-PSD-eigensolve budgets on identical integer starts."""
import subprocess,sys,json,hashlib,time,argparse
from pathlib import Path
here=Path(__file__).resolve().parent;source='candidates/refit100_legacy_exact_split_B479.npz'
p=argparse.ArgumentParser();p.add_argument('--prefix',default='');a=p.parse_args()
record=dict(input=source,source_sha256=hashlib.sha256((here/source).read_bytes()).hexdigest(),search_sha256=hashlib.sha256((here/'search.py').read_bytes()).hexdigest(),block_eigensolve_budget=650,cases=[])
for mode in ['single','coupled','thermal']:
    out=f'candidates/{a.prefix}matched_{mode}_B479.npz'
    cmd=[sys.executable,'search.py','--init',source,'--out',out,'--mode',mode,'--steps','1500','--seconds','3600','--seed','101','--refit-budget','650','--screen32']
    if mode=='thermal':cmd+=['--temperature','.000003','--subtree']
    start=time.time()
    with (here/f'results/{a.prefix}matched_{mode}_B479.log').open('w') as f:status=subprocess.run(cmd,cwd=here,stdout=f,stderr=subprocess.STDOUT).returncode
    record['cases'].append(dict(mode=mode,command=cmd,status=status,seconds=time.time()-start))
    (here/f'results/{a.prefix}matched_manifest.json').write_text(json.dumps(record,indent=2));print(json.dumps(record['cases'][-1]),flush=True)
