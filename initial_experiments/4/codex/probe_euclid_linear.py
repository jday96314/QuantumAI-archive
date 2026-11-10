"""Inspect the Euclidean cofactor transform as an in-place Sylvester factorization.

Read-only use of the attributed upstream classical model. No decoder claims.
"""
import importlib.util,inspect,json,random
from pathlib import Path
import numpy as np
from probe_wronskian import p,mul,model

path=Path('../../../initial_experiments/3/vendor/uncompute/tasks/dqi-reed-solomon/circuits/discovered/gen_constants.py')
spec=importlib.util.spec_from_file_location('classical_euclid',path);g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
src=inspect.getsource(g.backward)
src=src.replace('def backward(f, g, pre, m, log, hook=None, win=None):','def payload_R(f, g, pre, m, log, aa, bb, hook=None, win=None):')
src=src.replace('a, b = [0] * NH, [0] * NH\n    a[(m + 1) % NH] = ginv(f[17])','a, b = list(aa), list(bb)')
namespace=dict(g.__dict__);exec(src,namespace);payload_R=namespace['payload_R']

def trial(seed=81,L=35):
 rng=random.Random(seed);sig=[1]+[0]*35
 for pos in rng.sample(range(4095),L):
  root=p.EXP[pos]
  for k in reversed(range(1,36)):sig[k]^=mul(sig[k-1],root)
 E=[g.gsqrt(sig[2*i]) for i in range(18)];O=[g.gsqrt(sig[2*i+1]) for i in range(18)]
 f,h,eps,idx,pre,log=g.forward(E,O)
 K=np.zeros((35,36),dtype=np.int64);R=np.zeros((36,36),dtype=np.int64)
 for column in range(36):
  aa=[0]*18;bb=[0]*18
  if column<18:aa[column]=1
  else:bb[column-18]=1
  ee,oo,a,b=payload_R(f,h,pre,0,log,aa,bb)
  assert ee==E and oo==O
  R[:,column]=a+b
  for i in range(18):
   for j in range(18):K[i+j,column]^=mul(E[i],a[j])^mul(O[i],b[j])
 nz=np.count_nonzero(K,axis=0).tolist()
 print(json.dumps(dict(seed=seed,L=L,pre=pre,eps=eps,nonzero_per_column=nz,total_nonzero=sum(nz),K=K.tolist())),flush=True)
 return E,O,R,K

if __name__=='__main__':
 E,O,R,K=trial();np.savez('results/euclid_linear_probe.npz',E=E,O=O,R=R,K=K)
