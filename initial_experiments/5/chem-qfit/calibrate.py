import json,hashlib,shutil
import numpy as np
from core import *
from dfthcsos.verify import load_inputs,evaluate_words
from dfthcsos.scoreinfra.dfthc_math import reconstruct_x,build_h1_prime

z,u,w,tree,words=load(INC/'submission.npz');m=Metric(tree);ut=m.tensor(u);wt=m.tensor(w)
inputs=load_inputs(INC/'submission.npz');allwords=z['angle_words']
metrics=evaluate_words(inputs,allwords,tree,8)
decoded=reconstruct_u_vectors(allwords/128,tree);encoded=encode(decoded,tree,8,True)
roundtrip=evaluate_words(inputs,encoded,tree,8)
rec=dict(source=str(INC),metrics=metrics,roundtrip=roundtrip,
         word_mismatches=int((allwords!=encoded).sum()),
         tensor_metric=float(m.loss(ut,wt).sqrt()),eigen_min=float(m.ev.min()),
         projector_max_error=float(np.max(np.abs(np.einsum('ap,aq->apq',decoded,decoded)-np.einsum('ap,aq->apq',reconstruct_u_vectors(encoded/128,tree),reconstruct_u_vectors(encoded/128,tree))))),
         source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in INC.iterdir() if p.is_file()},
         controls=[])
for ridge in [1e-5,1e-8,0]:
    ww,err=m.project(ut,wt,1,ridge)
    # Also reproduce the unconditional ridge projection, rather than monotone guard.
    P=m.feats(ut)[0];S=P@P.T;L=torch.linalg.cholesky(S+ridge*torch.eye(len(P),**m.d));H=P@m.G@P.T
    H=torch.linalg.solve_triangular(L,H,upper=False);H=torch.linalg.solve_triangular(L,H.T,upper=False).T
    ev,V=torch.linalg.eigh((H+H.T)/2);wr=torch.linalg.solve_triangular(L.T,V[:,-350:]*ev[-350:].clamp(min=0).sqrt(),upper=True)[None]
    rec['controls'].append(dict(ridge=ridge,guarded_error=err,unconditional_error=float(m.loss(ut,wr).sqrt())))
for tag in ['L_r1b719c350','snap_B479','snap_B447','L_r2b351c175','qd_r1cal']:
    zz,uu,ww,tt,aw=load(ROOT/f'initial_experiments/4/claude/chem/{tag}.npz')
    mm=Metric(tt);uu=mm.tensor(uu);ww=mm.tensor(ww)
    entry=dict(tag=tag,shape=list(ww.shape),continuous_or_saved_error=float(mm.loss(uu,ww).sqrt()))
    if aw is None:aw=encode(uu.cpu().numpy().reshape(-1,76),tt,8,True)
    uq=mm.decode(mm.tensor(aw).reshape(*uu.shape[:-1],75));entry['direct_round_error']=float(mm.loss(uq,ww).sqrt())
    for ridge in [1e-5,0]:
        wwq,err=mm.project(uq,ww,4,ridge);entry[f'refit_ridge_{ridge}']=err
        if ridge==0:save(HERE/f'candidates/{tag}_rounded.npz',uq,wwq,tt,aw,err2=err)
    rec['controls'].append(entry)
    print(json.dumps(entry),flush=True)
out=HERE/'results/calibration.json';out.write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec),flush=True)
