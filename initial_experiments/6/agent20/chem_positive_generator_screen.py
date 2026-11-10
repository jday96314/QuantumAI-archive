"""Necessary positive-map witnesses for a nonnegative projector factorization.

If every generator Xc is PSD (column sign is irrelevant), then G=Σ|Xc><Xc|
maps every PSD test matrix Y to Σ Tr(Xc Y) Xc, which is PSD. For unit
v, a negative eigenvalue of unpack(G feat(v)) lower-bounds the Frobenius
error of every such factorization. This uses the fresh fixed native shift.
"""
import json,time,hashlib
from pathlib import Path
import numpy as np
import torch
from audit_core import Metric
from chem_rank_screen import BANK
HERE=Path(__file__).resolve().parent
if __name__=='__main__':
    torch.set_num_threads(1);started=time.time();z=dict(np.load(BANK));m=Metric(z['tree_edges'],device='cpu')
    rng=np.random.default_rng(2026102518);v=np.concatenate([np.eye(76),rng.normal(size=(256,76))]);v/=np.linalg.norm(v,axis=1,keepdims=True)
    f=m.feats(m.tensor(v));Y=m.unpack(f@m.G);ev,U=torch.linalg.eigh(Y);mins=ev[:,0].numpy();j=int(mins.argmin())
    bound=max(0.,-float(mins[j]));u=U[j,:,0].numpy();atom=m.feats(m.tensor(u[None]))[0]
    value=float(atom@(m.G@f[j]));assert abs(value-mins[j])<1e-10
    diagonal=np.flatnonzero(m.ii.numpy()==m.jj.numpy());off=np.flatnonzero(m.ii.numpy()!=m.jj.numpy())
    # G(offdiag,offdiag)/2 <= G(ii,jj) under PSD generators.
    G=m.G.numpy();idx={(int(i),int(k)):a for a,(i,k) in enumerate(zip(m.ii,m.jj))};viol=[]
    for a in off:
        i,k=int(m.ii[a]),int(m.jj[a]);viol.append(.5*G[a,a]-G[idx[(i,i)],idx[(k,k)]])
    row=dict(status='NECESSARY_PSD_GENERATOR_WITNESS_SCREEN',tests=len(v),minimum_output_eigenvalue=float(mins.min()),
        universal_error_lower_bound_for_PSD_generators=bound,witness_value=value,required_error=.02,
        witness_test_index=j,coordinate_tests_minimum=float(mins[:76].min()),max_pair_inequality_violation=float(max(viol)),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seconds=time.time()-started,
        scope='lower bound applies to all sums of PSD generator outer products under the fixed native shift; does not apply to unrestricted signed generators')
    np.savez(HERE/'parameters/positive_generator_witness.npz',input_v=v[j],output_v=u)
    (HERE/'results/positive_generator_screen.json').write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(row),flush=True)
