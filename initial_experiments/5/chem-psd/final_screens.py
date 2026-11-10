"""Wait for terminal geometry checkpoints, then evaluate decoded native physics.

Independent numerical subprocesses are computation, not delegated agents.
Only passing word screens build ANGLE families and invoke full verification.
"""
import concurrent.futures,json,os,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]

def wait_result(tag):
    path=HERE/f'results/{tag}.json'
    while not path.exists():time.sleep(5)
    return json.loads(path.read_text())

def screen(input_path,tag,C,steps=600,jacobi=0):
    env=dict(os.environ,PSD_DEVICE='cpu')
    cmd=[sys.executable,str(HERE/'physical.py'),'--input',str(input_path),'--tag',tag,'--rank',str(C),'--steps',str(steps),'--jacobi',str(jacobi),'--build']
    with (HERE/f'results/{tag}.log').open('w') as log:subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    rec=json.loads((HERE/f'results/{tag}_physical.json').read_text())
    if rec.get('native_verification_pending'):
        cmd=[sys.executable,str(ROOT/'from_host/verifiers/dfthcsos/verify.py'),rec['native_verification_pending'],'--salt','round5-chem-psd-'+tag]
        with (HERE/f'results/{tag}_native.json').open('w') as log:subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    return rec

def geometry(tag,C,steps):
    wait_result(tag)
    return screen(HERE/f'parameters/{tag}.npz','screen_'+tag,C,steps)

def integer():
    tags=['integer_claude_r2b351','integer_newton_r2b351','integer_projected_r2b351']
    data=[(wait_result(tag),tag) for tag in tags]
    best,tag=min(data,key=lambda x:x[0]['final']['direct_error_2'])
    rec=screen(HERE/f'parameters/{tag}.npz','screen_integer_best_c200',200,800,8)
    more=[]
    if best['final']['direct_error_2']<=.02 and not rec['passes_physical_gates']:
        for C in [208,224,256]:more.append(screen(HERE/f'parameters/{tag}.npz',f'screen_integer_best_c{C}',C,1200,8))
    return dict(selected_geometry=tag,physical=rec,expanded_rank_screens=more)

def main():
    jobs=[('trust_r4b224_cpu',128,600),('trust_prior_r2b351_cpu',192,600),('trust_eigen_r2b351',192,600),('trust_graph_r4b224',128,400)]
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        futures=[ex.submit(geometry,*j) for j in jobs]+[ex.submit(integer)]
        for f in concurrent.futures.as_completed(futures):
            results.append(f.result());(HERE/'results/final_screens.json').write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':main()
