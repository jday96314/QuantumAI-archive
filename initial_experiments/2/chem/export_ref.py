"""Dump the reference target matrix G (5776x5776) and the inputs the trainer needs."""
import numpy as np,sys
from dfthcsos.verify import load_inputs
i=load_inputs(sys.argv[1]);N=76
tei=np.asarray(i.ref['tei'],dtype=np.float64);S=i.S
G=tei.copy()
d=np.arange(N)
G[:,:,d,d]+=0.5*S[:,:,None];G[d,d,:,:]+=0.5*S[None,:,:]
np.savez(sys.argv[2],G=G.reshape(N*N,N*N).astype(np.float32),S=S,beta_1=i.beta_1,h1_thc=i.ref['h1_thc'],ecore=i.ref['ecore'],electrons=i.ref['electrons'],
         U=i.params['U'],w=i.w,W=i.W)
print('ok',np.linalg.norm(G))
