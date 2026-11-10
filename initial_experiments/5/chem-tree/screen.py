"""Matched tree-family screening with exact words and rank-C physical refits."""
import argparse
from core import *
def main():
    p=argparse.ArgumentParser();p.add_argument('--controls',nargs='+',default=['incumbent','b639_float','b351_float','b479_exact'])
    p.add_argument('--cycles',type=int,default=12);p.add_argument('--families',nargs='*');a=p.parse_args()
    H=target();manifest=[]
    for name in a.controls:
        path=SOURCES[name];z,u,w,inc,oldwords=load(path)
        source=dict(name=name,path=str(path),sha256=sha(path),shape=list(w.shape),continuous_error=physical_error(u,w,H))
        manifest.append(source);report(HERE/'results/sources.json',manifest)
        trees={**families(u,w,inc),**physical_families(H)};recordpath=HERE/f'results/screen_{name}.json'
        records=json.loads(recordpath.read_text()) if recordpath.exists() else []
        for family,tree in trees.items():
            if a.families and family not in a.families:continue
            start=time.time();words,v=quantize(u,tree);words=signed_words(words,tree);v=decode(words*np.pi/128,tree)
            score=surrogate(u,v,w,H)
            m=Refit(v,H);wr,err,hist=m.fit(w,a.cycles)
            out=HERE/f'candidates/{name}__{family}.npz'
            save(out,v,wr,tree,words,source=str(path),source_sha256=sha(path))
            rec=dict(control=name,family=family,shape=list(w.shape),**summary(tree),**score,
                     error_refit=err,refit_history=hist,cycles=a.cycles,seconds=time.time()-start,checkpoint=str(out))
            records=[r for r in records if r['family']!=family]+[rec];report(recordpath,records,False)
            print(json.dumps(rec),flush=True)
if __name__=='__main__':main()
