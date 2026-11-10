"""Search common rotation trees for low weighted quantization distortion."""
import argparse,heapq,json,time
from pathlib import Path
import numpy as np
from compress_chem import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--brot',type=int,default=8);p.add_argument('--steps',type=int,default=25000)
a=p.parse_args();root=Path(__file__).resolve().parent;z=dict(np.load(a.fit));u=z['U'][0];B,n=u.shape
base=z.get('tree_edges',np.load(root/'candidates/chem35c45_rot.npz')['tree_edges'])
K=z['K'] if 'K' in z else z['w'][0]@z['w'][0].T
S=(u@u.T)**2;weight=np.maximum(np.sum((K@S)*K,axis=1),1e-15);weight/=weight.mean()
alpha=np.pi/(1<<(a.brot-1));rng=np.random.default_rng(20261023);start=time.time()
def representative(node):return node if isinstance(node,int) else representative(node[0])
def edges(node):
    if isinstance(node,int):return []
    l,r=node;return [(representative(l),representative(r))]+edges(l)+edges(r)
def joined(l,r):
    if representative(r)==0:l,r=r,l
    return l,r
def huffman():
    mass=weight@(u*u);heap=[(float(v),i,i) for i,v in enumerate(mass)];heapq.heapify(heap);serial=n
    while len(heap)>1:
        w,_,l=heapq.heappop(heap);v,_,r=heapq.heappop(heap)
        heapq.heappush(heap,(w+v,serial,joined(l,r)));serial+=1
    return np.array(edges(heap[0][2]))
def greedy_quant():
    nodes=list(range(n));profiles=(u*u).T.copy()
    while len(nodes)>1:
        th=np.arctan2(np.sqrt(profiles[None,:,:]),np.sqrt(profiles[:,None,:]))
        dt=th-np.rint(th/alpha)*alpha
        cost=((profiles[:,None,:]+profiles[None,:,:])*dt*dt)@weight
        cost[np.tril_indices(len(nodes))]=np.inf;i,j=np.unravel_index(cost.argmin(),cost.shape)
        nd=joined(nodes[i],nodes[j]);pr=profiles[i]+profiles[j]
        keep=[k for k in range(len(nodes)) if k not in [i,j]]
        nodes=[nodes[k] for k in keep]+[nd];profiles=np.r_[profiles[keep],pr[None,:]]
    return np.array(edges(nodes[0]))
def evaluate(tree,detail=False):
    words=encode(u,tree,a.brot);uq=reconstruct_u_vectors(words/(1<<(a.brot-1)),tree)
    value=float(weight@(2*(1-np.sum(u*uq,axis=1)**2)))
    return (value,words,uq) if detail else value
trees=[('baseline',base),('huffman',huffman()),('greedy_quant',greedy_quant())]
for f in sorted((root/'parameters').glob('global500_tree_q8_*_tree.npy')):trees.append((f.stem,np.load(f)))
scores=[(evaluate(t),name,t) for name,t in trees];scores.sort(key=lambda v:v[0])
print(json.dumps(dict(initial=[(name,v) for v,name,_ in scores],seconds=time.time()-start)),flush=True)
best,_,tree=scores[0];tree=tree.copy();cur=best;besttree=tree.copy();accept=0
for step in range(a.steps):
    i,j=rng.choice(np.arange(1,n),size=2,replace=False)
    proposal=tree.copy();proposal[tree==i]=j;proposal[tree==j]=i
    value=evaluate(proposal);temp=scores[0][0]*.002*(.001+.999*(1-step/a.steps)**3)
    if value<cur or rng.random()<np.exp(min(0,(cur-value)/temp)):
        tree,cur=proposal,value;accept+=1
    if value<best:best,besttree=value,proposal.copy()
    if step%500==0 or step==a.steps-1:
        val,words,uq=evaluate(besttree,True)
        np.savez(a.out,U=u[None,:,:],tree_edges=besttree,err_proxy=val)
        np.savez(str(Path(a.out).with_suffix(''))+'_quantized.npz',U=uq[None,:,:],tree_edges=besttree,angle_words=words,brot=a.brot,err_proxy=val)
        print(json.dumps(dict(step=step,best=best,current=cur,accepted=accept,seconds=time.time()-start)),flush=True)
