"""Exact canonical-angle relaxation of two fixed block dictionary spaces.

For Q0^T Q1=U diag(c) V^T, arbitrary symmetric block coefficients have
independent two-by-two Gram systems with off-diagonal c_i*c_j. Solve these
systems (pseudoinverse on the exact common space) before imposing PSD/rank-C.
This is a relaxation of the actual model, not an achievable physical fit.
"""
import argparse,json,hashlib
from core import *
p=argparse.ArgumentParser();p.add_argument('fits',nargs='+');p.add_argument('--cpu',action='store_true');a=p.parse_args();rows=[]
for path in a.fits:
    z,u,w,tree,words=load(path);m=Metric(tree)
    if a.cpu:
        m.d['device']='cpu';m.G=m.G.cpu();m.fac=m.fac.cpu();m.ii=m.ii.cpu();m.jj=m.jj.cpu()
        m.stages=[tuple(x.cpu() for x in stage) for stage in m.stages]
    u=m.decode(m.tensor(words.reshape(*u.shape[:2],75)));P=m.feats(u);R,B,C=w.shape
    snapshot=HERE/f'candidates/capacity_snapshot_{Path(path).stem}.npz'
    save(snapshot,u,w,tree,words)
    with torch.no_grad():
        qs=[torch.linalg.qr(pp.T).Q for pp in P]
        if R==2:
            left,s,right=torch.linalg.svd(qs[0].T@qs[1]);E=qs[0]@left;F=qs[1]@right.T
            cc=s[:,None]*s[None,:];den=1-cc;common=den.abs()<1e-11
            def proj(G):
                h0=E.T@G@E;h1=F.T@G@F
                plus=(h0+h1)/(1+cc);minus=torch.where(common,torch.zeros_like(h0),(h0-h1)/den)
                aa=(plus+minus)/2;bb=(plus-minus)/2
                return E@aa@E.T+F@bb@F.T
            residual=m.G-proj(m.G)
            for _ in range(2):residual-=proj(residual)
            orth=max(float((q.T@residual@q).norm()) for q in qs)
            bound=float(residual.norm());nonzero=den[~common]
            row=dict(fit=path,shape=[R,B,C],two_block_symmetric_linear_relaxation=bound,normal_equation_residual=orth,common_principal_directions=int((s>1-1e-11).sum()),minimum_nonzero_Gram_eigenvalue=float(nonzero.abs().min()))
        else:
            q=qs[0];h=q.T@m.G@q;ev,V=torch.linalg.eigh(h);factor=q@(V[:,-C:]*ev[-C:].clamp(min=0).sqrt())
            row=dict(fit=path,shape=[R,B,C],fixed_dictionary_PSD_rankC_optimum=float((m.G-factor@factor.T).norm()))
        # SVD determines the actual union span, including exact duplicate rows.
        qq,sv,_=torch.linalg.svd(P.flatten(0,1).T,full_matrices=False);qq=qq[:,sv>sv[0]*1e-10]
        row['union_rank']=qq.shape[1];row['joint_span_relaxation']=float((m.G-qq@(qq.T@m.G@qq)@qq.T).norm())
        row['source_fit']=path;row['fit']=str(snapshot);row['snapshot_sha256']=hashlib.sha256(snapshot.read_bytes()).hexdigest()
        rows.append(row);print(json.dumps(row),flush=True)
(HERE/'results/capacity.json').write_text(json.dumps(rows,indent=2))
