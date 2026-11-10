"""Budgeted representation and joint exact-word/shared-tree searches."""
import argparse
from core import *

def main():
    p=argparse.ArgumentParser();p.add_argument('--control',default='b351_float');p.add_argument('--init')
    p.add_argument('--tag',required=True);p.add_argument('--steps',type=int,default=180)
    p.add_argument('--cycles',type=int,default=8);p.add_argument('--seed',type=int,default=179)
    p.add_argument('--mode',choices=['representation','random','joint'],default='representation')
    p.add_argument('--temperature',type=float,default=0.);a=p.parse_args()
    rng=np.random.default_rng(a.seed);H=target();path=Path(a.init) if a.init else SOURCES[a.control]
    z,physical,w0,inc,oldwords=load(path);R,B,C=w0.shape;start=time.time();history=[];seen=set();population=[]
    if a.mode=='joint':
        words=oldwords if oldwords is not None else quantize(physical,inc)[0].reshape(-1,75)
        words=words.reshape(R,B,75);u=decode(words*np.pi/128,inc);w,err,_=Refit(u,H).fit(w0,100)
        population=[dict(tree=inc,u=u,w=w,words=words,error=err)]
        initial_states=population.copy()
    else:
        for name,tr in families(physical,w0,inc,a.seed).items():
            sig=binary_signature(tr)
            if sig in seen:continue
            seen.add(sig);words,u=quantize(physical,tr);w,err,_=Refit(u,H).fit(w0,a.cycles)
            population.append(dict(tree=tr,u=u,w=w,words=words,error=err))
        population.sort(key=lambda x:x['error']);initial_states=population.copy();population=population[:8]
    best=population[0].copy();out=HERE/f'candidates/{a.tag}.npz'
    archive=HERE/'trajectories'/f'{a.tag}_dictionaries';archive.mkdir(exist_ok=True)
    def keep_dictionary(tr,words,err):
        file=archive/(fingerprint(tr)+'.npz')
        if file.exists() and float(np.load(file)['refit_error'])<=err:return
        np.savez_compressed(file,tree_edges=tr.astype(np.uint8),angle_words=words.astype(np.uint8),
                            brot=8,shape=np.array(w0.shape),refit_error=err)
    for s in initial_states:keep_dictionary(s['tree'],s['words'],s['error'])
    save(out,best['u'],best['w'],best['tree'],best['words'],err2=best['error'])
    current=best.copy();accepted=0;uphill=0;refits=0;best_origin='initial'
    for step in range(a.steps):
        mode=['orbital','siblings','reparent'][step%3]
        if a.mode=='random':tr=ordered(np.r_[0,rng.permutation(np.arange(1,N))]);source=physical;warm=w0
        elif a.mode=='representation':
            parent=population[int(rng.integers(min(4,len(population))))]
            tr=mutate(parent['tree'],rng,mode);source=physical;warm=w0
        else:
            # Include true matched fixed-word refitting before evaluating any
            # geometry change; weight-only progress cannot count as a move.
            current['w'],current['error'],_=Refit(current['u'],H).fit(current['w'],a.cycles)
            if current['error']<best['error']:
                best=current.copy();best_origin='fixed_word_refit'
                save(out,best['u'],best['w'],best['tree'],best['words'],err2=best['error'])
            tr=mutate(current['tree'],rng,mode);source=current['u'];warm=current['w']
        sig=binary_signature(tr)
        words,u=quantize(source,tr);words=signed_words(words,tr);u=decode(words*np.pi/128,tr)
        same=bool(np.max(np.abs(np.abs(np.sum(source*u,axis=-1))-1))<1e-13)
        if a.mode=='joint' and same:
            rec=dict(step=step,mode=mode,equivalent=True,current=current['error'],best=best['error'],refits=refits,seconds=time.time()-start)
        else:
            wr,err,hist=Refit(u,H).fit(warm,a.cycles);refits+=a.cycles*R
            keep_dictionary(tr,words,err)
            state=dict(tree=tr,u=u,w=wr,words=words,error=err)
            accept=False
            if a.mode=='joint':
                diff=err*err-current['error']**2
                temp=a.temperature*(.1+.9*(1-(step%40)/40))
                accept=bool(diff<0 or (temp>0 and rng.random()<np.exp(min(0,-diff/temp))))
                if accept:
                    uphill+=int(err>current['error']);current=state;accepted+=1
                if step%40==39:current=best.copy()
            else:
                if sig not in {binary_signature(x['tree']) for x in population}:
                    population.append(state);population.sort(key=lambda x:x['error']);population=population[:8]
            if err<best['error']:
                best=state.copy();best_origin='topology_proposal';save(out,u,wr,tr,words,err2=err)
                save(HERE/f'trajectories/{a.tag}_{step}.npz',u,wr,tr,words,err2=err)
            rec=dict(step=step,mode=mode,equivalent=same,proposal=err,best=best['error'],
                     current=current['error'],accepted=accept,uphill=uphill,refits=refits,best_origin=best_origin,
                     binary_signature=sig,seconds=time.time()-start,rng_state_after=rng.bit_generator.state)
        history.append(rec);report(HERE/f'results/{a.tag}.json',history,False)
        if step%10==0 or rec.get('accepted') or rec.get('proposal',1)<=best['error']:print(json.dumps(rec),flush=True)
    # Long exact-word final audit; preserve original short-budget search state.
    wr,err,hist=Refit(best['u'],H).fit(best['w'],300)
    save(HERE/f'candidates/{a.tag}_polished.npz',best['u'],wr,best['tree'],best['words'],err2=err)
    report(HERE/f'results/{a.tag}_final.json',dict(best_budget_error=best['error'],polished_error=err,
           cycles=a.cycles,steps=a.steps,source=str(path),mode=a.mode,seed=a.seed,refits=refits,
           polish_history=hist,seconds=time.time()-start))

if __name__=='__main__':main()
