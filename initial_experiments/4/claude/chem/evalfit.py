import sys,numpy as np,torch
sp=np.load('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters/reference_spectrum.npz')
ev,V=sp['eigenvalues'],sp['eigenvectors'];d=dict(device='cuda',dtype=torch.float64)
G=torch.tensor((V*np.maximum(ev,0))@V.T,**d);ii=torch.tensor(sp['ii'],device='cuda');jj=torch.tensor(sp['jj'],device='cuda');fac=torch.tensor(sp['fac'],**d)
for f in sys.argv[1:]:
    z=np.load(f);u=torch.tensor(z['U'],**d);w=torch.tensor(z['w'],**d)
    if u.dim()==2:u=u[None];w=w[None]
    un=u/u.norm(dim=-1,keepdim=True);F=un[...,ii]*un[...,jj]*fac
    X=torch.bmm(w.transpose(1,2),F).reshape(-1,len(ii));print(f,float((X.T@X-G).square().sum().sqrt()))
