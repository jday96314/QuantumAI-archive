"""Apply explicit finite research budgets; final PSD certification still runs.

Only targets PIDs whose command line contains this directory's geometry
script and the expected tag. The solver's SIGTERM handler finishes its saved
checkpoint and independent final convex solve.
"""
import json,os,signal,time
from pathlib import Path

HERE=Path(__file__).resolve().parent
tags=['trust_r4b224_cpu','trust_prior_r2b351_cpu','trust_eigen_r2b351','trust_graph_r4b224']
start=time.time();pending=set(tags);records=[]
while pending:
    for tag in list(pending):
        result=HERE/f'results/{tag}.json';progress=HERE/f'results/{tag}_progress.json'
        if result.exists():pending.remove(tag);continue
        try:data=json.loads(progress.read_text())
        except (FileNotFoundError,json.JSONDecodeError):continue
        n=len(data['history'])
        if n<100 and time.time()-start<1200:continue
        for proc in Path('/proc').iterdir():
            if not proc.name.isdigit():continue
            try:args=(proc/'cmdline').read_bytes().split(b'\0')
            except (OSError,ProcessLookupError):continue
            text=[s.decode(errors='replace') for s in args]
            if len(text)>2 and text[0].endswith('python') and any(s.endswith('chem-psd/geometry_trust.py') for s in text) and tag in text:
                os.kill(int(proc.name),signal.SIGTERM)
                records.append(dict(tag=tag,accepted_steps_at_signal=n,reason='100 accepted geometry steps or 1200-second continuation budget',pid=int(proc.name)))
                pending.remove(tag);break
        (HERE/'results/geometry_budgets.json').write_text(json.dumps(records,indent=2)+'\n')
    if pending:time.sleep(10)
