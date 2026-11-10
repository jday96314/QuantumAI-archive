"""Exact native pair metric and physical quantities; no verifier changes.

Reuses chem-qfit/core.py (exact block projections, gauge alignment), physical.py
(AGY offset derivatives), and round-4 Codex signed encoder. Host DFTHC math is
Apache-2.0, from from_host/verifiers/LICENSE.txt. All imports are read-only.
"""
from pathlib import Path
import json,os
import numpy as np
import torch
from core import Metric as PriorMetric, load, save, encode, reconstruct_u_vectors
from dfthcsos.verify import load_inputs, evaluate_words, passes_accuracy
from dfthcsos.scoreinfra.dfthc_math import build_h1_prime, reconstruct_x

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
INC=ROOT/'initial_experiments/4/codex/candidates/norm_cycle2_signed'
torch.set_num_threads(int(os.environ.get('AUDIT_THREADS','2')))

class Metric(PriorMetric):
    def __init__(self,tree,compare_cached=False,device=None):
        # Do not instantiate the prior cached-spectrum GPU target only to
        # replace it. This avoids an unnecessary 66 MiB device allocation.
        self.align_columns=True
        self.d=dict(device=device or os.environ.get('AUDIT_DEVICE','cuda'),dtype=torch.float64)
        sp=np.load(ROOT/'initial_experiments/3/parameters/reference_spectrum.npz')
        self.ii=torch.tensor(sp['ii'],device=self.d['device']);self.jj=torch.tensor(sp['jj'],device=self.d['device'])
        self.fac=self.tensor(sp['fac']);self.ev=sp['eigenvalues'];self.tree=tree;self.stages=[]
        for j,(p,c) in enumerate(tree):
            if not self.stages or any(p in (pp,cc) or c in (pp,cc) for _,pp,cc in self.stages[-1]):self.stages.append([])
            self.stages[-1].append((j,int(p),int(c)))
        self.stages=[tuple(torch.tensor(x,device=self.d['device']) for x in zip(*stage)) for stage in self.stages]
        self.inp=load_inputs(INC/'submission.npz')
        ii=self.ii.cpu().numpy();jj=self.jj.cpu().numpy();fac=self.fac.cpu().numpy()
        # Construct target directly from the freshly checked native FCIDUMP and
        # float32-rounded fixed shift. Packed coordinates are orthonormal.
        G=self.inp.ref['tei'][ii[:,None],jj[:,None],ii[None,:],jj[None,:]].copy()
        diag=(ii==jj).astype(float);s=self.inp.S[ii,jj]
        G+=.5*(s[:,None]*diag[None,:]+diag[:,None]*s[None,:])
        G*=fac[:,None]*fac[None,:]
        self.cached_target_difference=float(np.linalg.norm(G-(sp['eigenvectors']*self.ev)@sp['eigenvectors'].T)) if compare_cached else None
        self.G=self.tensor(G)
        self.h0=self.tensor(build_h1_prime(self.inp.ref['h1_thc'],self.inp.S,self.inp.beta_1,
            self.inp.ref['electrons'],np.zeros_like(self.inp.W),reconstruct_x(self.inp.params['U'],self.inp.w)))
        self.constant=float(self.inp.energy_upper-self.inp.ref['ecore']+self.inp.ref['electrons']*self.inp.beta_1)

    def unpack(self,x):
        v=x/self.fac
        h=torch.zeros((*v.shape[:-1],76,76),**self.d)
        h[...,self.ii,self.jj]=v;h[...,self.jj,self.ii]=v
        return h

    def physical(self,u,w,off):
        u=u/u.norm(dim=-1,keepdim=True)
        coefficients=(w@off[...,None]).squeeze(-1)
        h=self.h0-torch.bmm(u.transpose(1,2),u*coefficients[...,None]).sum(0)
        eig=torch.linalg.eigvalsh(h)
        absolute=eig.abs().sum()
        gap=self.constant+absolute-eig.sum()+.5*off.square().sum()
        lam=absolute+.25*(w.abs().sum(1)+(off+w.sum(1)).abs()).square().sum()
        return gap,lam,gap*(2*lam-gap)

    @torch.no_grad()
    def project_gap(self,u,w,off,target=3.91515):
        X=(w.transpose(1,2)@self.feats(u)).flatten(0,1)
        for _ in range(60):
            eig,V=torch.linalg.eigh(self.h0-self.unpack(off.flatten()@X))
            gap=self.constant+eig.abs().sum()-eig.sum()+.5*off.square().sum()
            delta=target-float(gap)
            if abs(delta)<1e-10:break
            N=V[:,eig<0]@V[:,eig<0].T
            dual=N[self.ii,self.jj]*self.fac
            grad=off.flatten()+2*X@dual
            update=(delta*grad/grad.square().sum().clamp(min=1e-20)).reshape(off.shape)
            update*=min(1.,2./max(float(update.norm()),1e-30))
            old=off.clone();changed=False
            for scale in [1.,.5,.25,.125,.0625,.015625,.00390625]:
                off.copy_(old+scale*update)
                newgap=self.physical(u,w,off)[0]
                if abs(float(newgap)-target)<abs(delta):changed=True;break
            if not changed:off.copy_(old);break
        return float(self.physical(u,w,off)[0])

def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def checkpoint(m,u,w,off,tree,cycles=2):
    """Float, direct rounded, exact refit metrics including normalization.

    Offsets start in the same aligned column gauge. Refit gap is projected to
    the floor margin, explicitly reported, without changing error_2.
    """
    with torch.no_grad():
        words=encode(u.cpu().numpy().reshape(-1,76),tree,8,True)
        uq=m.decode(m.tensor(words).reshape(*u.shape[:-1],75))
        wq,er=m.project(uq,w,cycles)
        oq=off.clone();m.project_gap(uq,wq,oq)
        def row(uu,ww,oo):
            g,l,s=m.physical(uu,ww,oo)
            domain=bool(float(g)>0 and float(2*l-g)>0)
            return dict(error_2=float(m.loss(uu,ww).sqrt()),Lambda=float(l),E_gap=float(g),
                lambda_eff=float(s.sqrt()) if domain else None,gap_in_domain=domain,
                weight_norm=float(ww.norm()),max_weight=float(ww.abs().max()))
        return dict(float=row(u,w,off),rounded=row(uq,w,off),refit=row(uq,wq,oq)),(uq,wq,oq,words)
