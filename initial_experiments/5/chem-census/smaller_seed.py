"""Continue measured one-row deletion from B511 to a fair B479 growth seed."""
from census_common import *
z,u,w,tree,words=load(HERE/'parameters/prune_B511.npz');m=Metric(tree)
u=m.tensor(u);w=m.tensor(w);off=m.tensor(z['W']);indices=z['indices']
with torch.no_grad():
    while u.shape[1]>479:
        P=m.feats(u)[0];S=P@P.T;iv=torch.linalg.inv(S);K=w[0]@w[0].T;d=iv.diag()
        score=2*(K@S@K).diag()/d-(K.diag()/d).square()
        ix=torch.arange(u.shape[1],device=u.device);ix=ix[ix!=score.argmin()]
        u=u[:,ix];w=w[:,ix];words=words[ix.cpu().numpy()];indices=indices[ix.cpu().numpy()]
        w,_=m.project(u,w,1)
    store('prune_B479',m,u,w,off,words,indices=indices)
    rec=quantities(m,u,w,off);write('smaller_seed.json',rec);print(json.dumps(rec))
