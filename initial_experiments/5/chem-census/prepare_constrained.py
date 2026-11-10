"""Transfer banked physical generators to the *repaired* B703 dictionary."""
from census_common import *
z,m,u,w,off,words=initial()
with torch.no_grad():
    X=w[0].T@m.feats(u)[0]
    zz,uu,ww,_,aw=load(HERE/'parameters/repair_B703_extended.npz')
    ut=m.tensor(uu);P=m.feats(ut)[0]
    wt=torch.linalg.lstsq(P.T,X.T).solution[None]
    store('continuation703_physical_seed',m,ut,wt,off,aw)
    rec=quantities(m,ut,wt,off)
    write('continuation703_physical_seed.json',rec);print(json.dumps(rec))
