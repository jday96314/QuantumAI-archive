import sys, numpy as np, torch
sys.path.insert(0,'/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3')
sys.path.insert(0,'/media/shared/ExpansionDrive2/QuantumAI/from_host/verifiers')
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
z=np.load(sys.argv[1]);tree=z['tree_edges'];words=z['angle_words'];brot=int(z['brot'])
R,B,C=z['w'].shape
sp=np.load('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/reference_spectrum.npz')
ev,V=sp['eigenvalues'],sp['eigenvectors'];d=dict(device='cuda',dtype=torch.float64)
G=torch.tensor((V*np.maximum(ev,0))@V.T,**d);ii=torch.tensor(sp['ii'],device='cuda');jj=torch.tensor(sp['jj'],device='cuda');fac=torch.tensor(sp['fac'],**d)
def proj(u):
    u=torch.tensor(u,**d);F=u[:,ii]*u[:,jj]*fac
    Q,Rm=torch.linalg.qr(F.T);Hp=Q.T@G@Q;e,v=torch.linalg.eigh((Hp+Hp.T)/2)
    idx=torch.argsort(e,descending=True)[:C];root=v[:,idx]*e[idx].clamp(min=0).sqrt()
    w=torch.linalg.solve(Rm,root);X=w.T@F;return float((X.T@X-G).square().sum().sqrt())
for bits in [brot,brot-1,brot-2]:
    sh=brot-bits;wq=(((words>>sh)+((words>>(sh-1))&1 if sh else 0))%(1<<bits)).astype(np.int64)
    u=reconstruct_u_vectors(wq/(1<<(bits-1)),tree)
    print(bits,proj(u),flush=True)
