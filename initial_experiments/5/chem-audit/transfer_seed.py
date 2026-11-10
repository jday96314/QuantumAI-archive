"""Project selected physical generators into the improved integer dictionary."""
from audit_core import *
m=Metric(np.load(INC/'submission.npz')['tree_edges']);base=dict(np.load(INC/'submission.npz'))
u=m.decode(m.tensor(base['angle_words'][76:]).reshape(1,719,75));P=m.feats(u)
z,un,wn,tree,words=load(ROOT/'initial_experiments/5/chem-qfit/candidates/discrete_control.npz')
un=m.decode(m.tensor(words).reshape(1,719,75));Pn=m.feats(un)
w=m.tensor(m.inp.w);off=m.tensor(m.inp.W)
with torch.no_grad():
    F=w.transpose(1,2)@P
    ww=torch.linalg.solve(Pn@Pn.transpose(1,2),Pn@F.transpose(1,2))
    m.project_gap(un,ww,off)
    save(HERE/'parameters/transferred_seed.npz',un,ww,tree,words,W=off.cpu().numpy())
    rec=checkpoint(m,un,ww,off,tree)[0]
    write_json(HERE/'results/transferred_seed.json',rec);print(json.dumps(rec),flush=True)
