"""Bound the expensive lookahead pilot, then independently certify its NPZ."""
import os,signal,time
from core import *
from solve import solve

tag='integer_projected_r2b351';path=HERE/f'results/{tag}.log'
while True:
    records=[]
    for line in path.read_text().splitlines():
        try:rec=json.loads(line)
        except json.JSONDecodeError:continue
        if 'step' in rec:records.append(rec)
    if records and records[-1]['step']>=10:break
    time.sleep(3)
pid=None
for proc in Path('/proc').iterdir():
    if not proc.name.isdigit():continue
    try:args=[s.decode(errors='replace') for s in (proc/'cmdline').read_bytes().split(b'\0')]
    except OSError:continue
    if args and args[0].endswith('python') and any(s.endswith('chem-psd/integer_trust.py') for s in args) and tag in args:
        pid=int(proc.name);os.kill(pid,signal.SIGINT);break
if pid:
    while Path(f'/proc/{pid}').exists():time.sleep(1)
z=dict(np.load(HERE/f'parameters/{tag}.npz'));m=Covariance(z['U'])
with torch.no_grad():
    S,h=solve(m,tensor(z['S']),2000,None,'fista',200,quiet=True)
    final=m.diagnostics(S,True);z.update(S=S.cpu().numpy(),w=m.weights(S).cpu().numpy())
    np.savez(HERE/f'parameters/{tag}.npz',**z)
    report(HERE/f'results/{tag}.json',dict(input=str(HERE/'parameters/integer_newton_seed.npz'),
        initial=dict(direct_error_2=records[0]['previous_error_2']),final=final,history=records,
        history_contains_logged_iterations_only=True,accepted_iterations_at_least=records[-1]['step']+1,
        termination='Finite projected-trial pilot budget; independent terminal PSD solve completed',final_history=h,
        device=D['device'],stopped_pid=pid))
