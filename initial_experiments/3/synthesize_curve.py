"""Search shared Davio expansions for the upper coefficient bits."""
import itertools,json
from pathlib import Path
import numpy as np

def network(table,outputs,order):
    n=len(order);x=np.arange(1<<n);ix=np.zeros_like(x)
    for k,j in enumerate(order):ix|=((x>>k)&1)<<j
    cache={};nodes=[]
    def visit(v,k):
        flip=int(v[0]);v=v^flip
        if not np.any(v):return ('c',flip)
        key=(k,np.packbits(v).tobytes())
        if key not in cache:
            half=len(v)//2;lo=v[:half];delta=lo^v[half:]
            aa=visit(lo,k-1);bb=visit(delta,k-1)
            if bb==('c',0):res=aa
            elif bb==('c',1):res=('x',aa,('v',order[k]))
            else:
                res=('x',aa,('a',len(nodes)));nodes.append((order[k],bb))
            cache[key]=res
        res=cache[key]
        return ('x',res,('c',1)) if flip else res
    roots=[visit(((table[ix]>>b)&1).astype(np.uint8),n-1) for b in outputs]
    return nodes,roots

if __name__=='__main__':
 root=Path(__file__).resolve().parent;T=np.load(root/'candidates/margin_curve8_t5.9e-08.npy')[:,2]
 assert np.array_equal(T[::2]>>4,T[1::2]>>4)
 best=1000
 for k,perm in enumerate(itertools.permutations(range(1,9))):
  # Bit zero is absent from the high coefficient words.
  nodes,roots=network(T,range(4,7),perm)
  if len(nodes)<best:
   best=len(nodes);rec=dict(ands=best,order=perm,nodes=nodes,roots=roots)
   (root/'parameters/curve_davio.json').write_text(json.dumps(rec));print(json.dumps(dict(iteration=k,ands=best,order=perm)),flush=True)
