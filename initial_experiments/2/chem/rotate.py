"""Minimize the two-body one-norm by per-rho orthogonal mixing of the C columns of (w, W).

Sum_c X_c (x) X_c, sum_c W_c X_c and sum W^2 are invariant, so error_2, h1', E_SOS and
the angle words are unchanged; only lambda_h2 = 1/4 sum (|w_c|_1 + |W_c + sum_b w_bc|)^2 moves.
"""
import sys,numpy as np,torch
src,dst=sys.argv[1],sys.argv[2];steps=int(sys.argv[3]) if len(sys.argv)>3 else 3000
z=dict(np.load(src));w=z['W'];Wu=z['H2_I']-w.sum(1);R,B,C=w.shape
dev='cuda';d=dict(dtype=torch.float64,device=dev)
wt=torch.tensor(w,**d);Wt=torch.tensor(Wu,**d)
A=torch.zeros(R,C,C,**d,requires_grad=True)
def lam(O):
    w2=torch.einsum('rbc,rcd->rbd',wt,O);W2=torch.einsum('rc,rcd->rd',Wt,O)
    return 0.25*((w2.abs().sum(1)+(W2+w2.sum(1)).abs())**2).sum()
opt=torch.optim.Adam([A],lr=2e-3)
best=(lam(torch.eye(C,**d).expand(R,C,C)).item(),None);print('start',best[0],flush=True)
for s in range(steps):
    opt.zero_grad();O=torch.linalg.matrix_exp(A-A.transpose(1,2));l=lam(O);l.backward();opt.step()
    if l.item()<best[0]:best=(l.item(),O.detach().clone())
    if s%500==0:print(s,l.item(),flush=True)
O=best[1]
w2=torch.einsum('rbc,rcd->rbd',wt,O).cpu().numpy();W2=torch.einsum('rc,rcd->rd',Wt,O).cpu().numpy()
z['W']=w2;z['H2_I']=W2+w2.sum(1);np.savez(dst,**z);print('final lambda_h2',best[0])
