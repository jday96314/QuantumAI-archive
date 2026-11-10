"""Eight-bit-spaced codebooks with global edge phase, legally emitted at 9 bits.

There are still 128 signed representatives per edge, just as in the incumbent
signed eight-bit representation. The extra low output bit is constant over
addresses. Host brot=9 is explicit and all additional costs must be charged.
No nonuniform semantics are silently assigned to eight-bit host words.
"""
import argparse
from core import *

def phase_encode(u,tree,phase,bits=9):
    scale=1<<(bits-1);stride=1<<(bits-8)
    v=np.array(u,copy=True);words=np.zeros((*u.shape[:-1],75),dtype=np.int64)
    for j in reversed(range(75)):
        p,c=tree[j];a=v[...,p].copy();b=v[...,c].copy();sign=np.where(a<0,-1.,1.)
        angle=np.arctan2(b*sign,a*sign)*scale/np.pi
        q=np.rint((angle-phase[j])/stride).astype(np.int64)*stride+phase[j]
        flip=q>=scale//2;q=np.where(flip,q-scale,q);sign=np.where(flip,-sign,sign)
        flip=q<-(scale//2);q=np.where(flip,q+scale,q);sign=np.where(flip,-sign,sign)
        words[...,j]=q%(2*scale);v[...,p]=sign*np.hypot(a,b);v[...,c]=0
    return words,decode(words*np.pi/scale,tree)

def main():
    p=argparse.ArgumentParser();p.add_argument('--control',default='b351_float');p.add_argument('--tree-init')
    p.add_argument('--tag',required=True);p.add_argument('--passes',type=int,default=2);p.add_argument('--bits',type=int,default=9)
    p.add_argument('--cycles',type=int,default=12);p.add_argument('--seed',type=int,default=87);a=p.parse_args()
    z,u,w,tree,_=load(SOURCES[a.control]);H=target();rng=np.random.default_rng(a.seed)
    if a.tree_init:tree=load(a.tree_init)[3]
    phase=np.zeros(75,dtype=np.int64);words,v=phase_encode(u,tree,phase,a.bits)
    wr,err,_=Refit(v,H).fit(w,a.cycles);start=time.time();history=[];initial=err;accepted=0;evals=1
    stride=1<<(a.bits-8);out=HERE/f'candidates/{a.tag}.npz'
    def checkpoint():np.savez(out,U=v,w=wr,tree_edges=tree,angle_words=words.reshape(-1,75),brot=a.bits,grid_phase=phase,err2=err)
    checkpoint()
    for repeat in range(a.passes):
        for j in rng.permutation(75):
            best=err;state=None
            for offset in range(stride):
                if offset==phase[j]:continue
                new=phase.copy();new[j]=offset;nwords,nu=phase_encode(u,tree,new,a.bits)
                nw,ne,hist=Refit(nu,H).fit(w,a.cycles);evals+=1
                if ne<best:best=ne;state=(new,nwords,nu,nw)
            changed=state is not None
            if changed:phase,words,v,wr=state;err=best;accepted+=1;checkpoint()
            rec=dict(pass_index=repeat,edge=int(j),error=err,changed=changed,accepted=accepted,
                     evaluations=evals,block_eigensolves=evals*a.cycles*w.shape[0],seconds=time.time()-start)
            history.append(rec);report(HERE/f'results/{a.tag}.json',history,False)
            if len(history)%15==0:print(json.dumps(rec),flush=True)
    wr,err,polish=Refit(v,H).fit(wr,300);checkpoint()
    report(HERE/f'results/{a.tag}_final.json',dict(initial=initial,final=err,passes=a.passes,
        bits=a.bits,phase=phase.tolist(),evaluations=evals,accepted=accepted,polish_history=polish,
        seconds=time.time()-start,semantics='host integer uniform grid, shared edge phases, charged brot'))

if __name__=='__main__':main()
