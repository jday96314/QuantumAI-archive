"""Replay frozen search decisions to archive every evaluated integer dictionary.

No numerical search or new score selection occurs. Recovered tree signatures
must match the saved trajectory at every step. Float64 best weights are saved
by the main search; these compact archives preserve every exact integer table.
"""
import argparse
from core import *

def main():
    p=argparse.ArgumentParser();p.add_argument('--control',default='b351_float');p.add_argument('--tag',required=True)
    p.add_argument('--mode',choices=['representation','random','joint'],default='representation')
    p.add_argument('--cycles',type=int,default=12);p.add_argument('--seed',type=int,default=179);a=p.parse_args()
    history=json.loads((HERE/f'results/{a.tag}.json').read_text());z,physical,w,inc,oldwords=load(SOURCES[a.control])
    rng=np.random.default_rng(a.seed);H=target();out=HERE/'trajectories'/f'{a.tag}_dictionaries';out.mkdir(exist_ok=True)
    def state(tr,u,words,error):return dict(tree=tr,u=u,words=words,error=error)
    def keep(s):
        file=out/(fingerprint(s['tree'])+'.npz')
        if file.exists() and float(np.load(file)['refit_error'])<=s['error']:return
        np.savez_compressed(file,tree_edges=s['tree'].astype(np.uint8),angle_words=s['words'].astype(np.uint8),
                            brot=8,shape=np.array(w.shape),refit_error=s['error'])
    population=[];seen=set()
    if a.mode=='joint':
        words=oldwords.reshape(*physical.shape[:-1],75);u=decode(words*np.pi/128,inc)
        _,err,_=Refit(u,H).fit(w,100);population=[state(inc,u,words,err)]
    else:
        for name,tr in families(physical,w,inc,a.seed).items():
            sig=binary_signature(tr)
            if sig in seen:continue
            seen.add(sig);words,u=quantize(physical,tr);_,err,_=Refit(u,H).fit(w,a.cycles)
            s=state(tr,u,words,err);population.append(s);keep(s)
        population.sort(key=lambda x:x['error']);population=population[:8]
    best=population[0].copy();current=best.copy();keep(best)
    previous_uphill=0
    for rec in history:
        step=rec['step'];mode=['orbital','siblings','reparent'][step%3]
        if a.mode=='random':tr=ordered(np.r_[0,rng.permutation(np.arange(1,N))]);source=physical
        elif a.mode=='representation':
            parent=population[int(rng.integers(min(4,len(population))))];tr=mutate(parent['tree'],rng,mode);source=physical
        else:tr=mutate(current['tree'],rng,mode);source=current['u']
        words,u=quantize(source,tr);words=signed_words(words,tr);u=decode(words*np.pi/128,tr)
        sig=binary_signature(tr)
        if 'binary_signature' in rec:assert sig==rec['binary_signature'],(step,sig,rec['binary_signature'])
        if 'proposal' in rec:
            s=state(tr,u,words,rec['proposal']);keep(s)
            if a.mode=='joint':
                if rec['best']<best['error']:
                    if rec['best_origin']=='topology_proposal':best=s.copy()
                    else:best=current.copy();best['error']=rec['best']
                if rec['accepted']:current=s
                # A rejected proposal necessarily used an annealing draw.
                # Accepted uphill draws are identified by the monotone
                # uphill counter; accepted descending proposals use none.
                if not rec['accepted'] or rec['uphill']>previous_uphill:rng.random()
                previous_uphill=rec['uphill']
                if step%40==39:current=best.copy()
            else:
                if sig not in {binary_signature(x['tree']) for x in population}:
                    population.append(s);population.sort(key=lambda x:x['error']);population=population[:8]
    report(HERE/f'results/{a.tag}_archive.json',dict(steps=len(history),exact_dictionaries=len(list(out.glob('*.npz'))),
           all_tree_signatures_match=True,directory=str(out),note='Exact words; main NPZs preserve float64 best weights'))
if __name__=='__main__':main()
