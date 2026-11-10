"""Choose a common Givens tree using vector co-occurrence and quantization loss."""
import argparse,json
from pathlib import Path
import numpy as np
from compress_chem import encode
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
a=p.parse_args();z=np.load(a.fit);u=z['U'][0];n=u.shape[1];rng=np.random.default_rng(441)
def edges(node):
 if isinstance(node,int):return []
 left,right=node
 return [(rep(left),rep(right))]+edges(left)+edges(right)
def rep(node):return node if isinstance(node,int) else rep(node[0])
def make(kind,jitter=0):
 nodes=list(range(n));profiles=(u*u).T.copy()
 while len(nodes)>1:
  m=len(nodes);pairs=[];unused=set(range(m))
  if kind=='overlap':score=np.minimum(profiles[:,None,:],profiles[None,:,:]).sum(2)
  elif kind=='cosine':score=np.sqrt(profiles)@np.sqrt(profiles).T
  elif kind=='random':score=rng.random((m,m))
  else:
   th=np.arctan2(np.sqrt(profiles[None,:,:]),np.sqrt(profiles[:,None,:]))
   dt=th-np.rint(th*128/np.pi)*np.pi/128
   score=-((profiles[:,None,:]+profiles[None,:,:])*dt*dt).sum(2)
  score+=rng.normal(size=(m,m))*jitter*max(1e-10,abs(score).max())
  new=[];pr=[]
  while len(unused)>1:
   ix=np.array(sorted(unused));sub=score[np.ix_(ix,ix)].copy();np.fill_diagonal(sub,-np.inf)
   q,r=np.unravel_index(sub.argmax(),sub.shape);i,j=int(ix[q]),int(ix[r]);unused.remove(i);unused.remove(j)
   l,rr=nodes[i],nodes[j]
   if rep(rr)==0:l,rr=rr,l
   new.append((l,rr));pr.append(profiles[i]+profiles[j])
  if unused:
   i=unused.pop();new.append(nodes[i]);pr.append(profiles[i])
  nodes,profiles=new,np.array(pr)
 return np.array(edges(nodes[0]),dtype=int)
trees=[('baseline',z['tree_edges'])]
for kind in ['overlap','cosine','quant','random']:
 for t in range(8):trees.append((kind+str(t),make(kind,0 if t==0 else .03)))
for bits in [8,9,10,11]:
 ranking=[]
 for name,tree in trees:
  words=encode(u,tree,bits);uq=reconstruct_u_vectors(words/(1<<(bits-1)),tree)
  error=float(np.square(u-uq).sum());ranking.append((error,name,tree,uq,words))
 ranking.sort(key=lambda x:x[0]);print(json.dumps(dict(brot=bits,ranking=[(x[1],x[0]) for x in ranking[:6]])),flush=True)
 for k,(_,name,tree,uq,words) in enumerate(ranking[:3]):
  out=a.out+f'_q{bits}_{k}'
  np.save(out+'_tree.npy',tree)
  np.savez(out+'.npz',U=uq[None,:,:],angle_words=words,tree_edges=tree,brot=bits)
