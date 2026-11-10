"""Certify unrestricted block-PSD relaxations of every screened decoded model."""
import concurrent.futures,json,time
from core import *
from solve import solve

def process(path):
    tag=path.name.removesuffix('_physical.json');out=HERE/f'results/{tag}_decoded_psd.json'
    if out.exists():return json.loads(out.read_text())
    z=np.load(HERE/f'parameters/{tag}_quantized.npz');m=Covariance(z['U']);S0=tensor(z['S'])
    with torch.no_grad():
        S,hist=solve(m,S0,1600,None,'fista',200,quiet=True)
        diag=m.diagnostics(S,True)
        np.savez(HERE/f'parameters/{tag}_decoded_psd.npz',U=z['U'],S=S.cpu().numpy(),w=m.weights(S).cpu().numpy(),tree_edges=z['tree_edges'],angle_words=z['angle_words'])
        rec=dict(physical_screen=str(path),diagnostics=diag,history=hist)
        report(out,rec)
    return rec

def main():
    seen=set();futures=[];summary=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        while True:
            for path in (HERE/'results').glob('screen_*_physical.json'):
                if path.name not in seen:seen.add(path.name);futures.append(pool.submit(process,path))
            finished=[f for f in futures if f.done()]
            for f in finished:summary.append(f.result());futures.remove(f)
            (HERE/'results/decoded_psd_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
            complete=HERE/'results/final_screens.json'
            if complete.exists() and len(json.loads(complete.read_text()))==5 and not futures:break
            time.sleep(5)
if __name__=='__main__':main()
