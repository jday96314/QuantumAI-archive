"""Enumerated local binary topology repair from exact word dictionaries."""
import argparse
from core import *

def main():
    p=argparse.ArgumentParser();p.add_argument('--control',default='b479_exact');p.add_argument('--init')
    p.add_argument('--tag',required=True);p.add_argument('--rounds',type=int,default=2)
    p.add_argument('--support',type=int,default=8);p.add_argument('--cycles',type=int,default=12);a=p.parse_args()
    path=Path(a.init) if a.init else SOURCES[a.control];z,u,w,tree,words=load(path);R,B,C=w.shape;H=target()
    w,err,_=Refit(u,H).fit(w,200);initial=err;start=time.time();history=[];total=0;accepted=0
    folder=HERE/'trajectories'/f'{a.tag}_dictionaries';folder.mkdir(exist_ok=True);out=HERE/f'candidates/{a.tag}.npz'
    save(out,u,w,tree,words,err2=err);top=[]
    for repeat in range(a.rounds):
        # Long fixed-word baseline prevents slow coefficient convergence from
        # masquerading as a topology improvement.
        w,err,_=Refit(u,H).fit(w,200);best=(err,u,w,tree,words)
        for i,tr in enumerate(nni_neighbors(tree,a.support)):
            nw,nu=quantize(u,tr);nw=signed_words(nw,tr);nu=decode(nw*np.pi/128,tr)
            ww,ne,hist=Refit(nu,H).fit(w,a.cycles);total+=1
            name=fingerprint(tr);np.savez_compressed(folder/(name+'.npz'),tree_edges=tr.astype(np.uint8),angle_words=nw.astype(np.uint8),brot=8,shape=np.array(w.shape),refit_error=ne)
            row=dict(round=repeat,index=i,error=ne,baseline=err,seconds=time.time()-start,signature=binary_signature(tr))
            history.append(row);report(HERE/f'results/{a.tag}.json',history,False)
            if ne<best[0]:best=(ne,nu,ww,tr,nw)
            if repeat==0:
                top.append((ne,name))
                top.sort();top=top[:3]
                if name in {x[1] for x in top}:save(HERE/f'candidates/{a.tag}_trial_{name}.npz',nu,ww,tr,nw,err2=ne)
            if i%15==0:print(json.dumps(row),flush=True)
        if best[0]>=err-1e-10:break
        err,u,w,tree,words=best;accepted+=1;save(out,u,w,tree,words,err2=err)
    # Keep all top candidates named independently of graph hashes for repair.
    for i,(value,name) in enumerate(top):
        data=dict(np.load(HERE/f'candidates/{a.tag}_trial_{name}.npz'));np.savez(HERE/f'candidates/{a.tag}_top{i}.npz',**data)
    wr,err,polish=Refit(u,H).fit(w,300);save(out,u,wr,tree,words,err2=err)
    report(HERE/f'results/{a.tag}_final.json',dict(initial=initial,final=err,accepted_topologies=accepted,
           evaluated_topologies=total,refit_cycles_per_proposal=a.cycles,top_candidates=top,
           support=a.support,seconds=time.time()-start,polish_history=polish))
if __name__=='__main__':main()
