"""Local outputs; read-only native chemistry and prior exact-word helpers.

Native DFTHC mathematics/costs are Apache-2.0 (from_host/verifiers/LICENSE.txt).
Exact pair projection and signed encoder reuse chem-qfit/core.py and round-4
Codex angle_gauge_probe.py. Fresh native target/physics reuse chem-audit.
"""
from pathlib import Path
import hashlib,json,time
import numpy as np
import torch
from audit_core import Metric,load,save,encode,reconstruct_u_vectors
from dfthcsos.verify import load_inputs,evaluate_words,passes_accuracy,heuristic_costs,resource_costs
from dfthcsos.scoreinfra.dfthc_math import reconstruct_x,build_h1_prime
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BANK=ROOT/'initial_experiments/5/chem-audit/candidates/banked'
BASE_COST=708156073.7876049
BASE_LAMBDA=48.5673374743
TARGET=BASE_COST/1.25
torch.set_num_threads(2)
def write(name,value):
    (HERE/'results'/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def initial():
    z=dict(np.load(BANK/'submission.npz'));tree=z['tree_edges'];m=Metric(tree)
    words=z['angle_words'][76:];u=m.decode(m.tensor(words)[None]);w=m.tensor(z['W']);off=m.tensor(z['H2_I']-z['W'].sum(1))
    return z,m,u,w,off,words
def store(tag,m,u,w,off,words,**kw):
    save(HERE/'parameters'/f'{tag}.npz',u,w,m.tree,words,W=off.detach().cpu().numpy(),**kw)
def quantities(m,u,w,off):
    with torch.no_grad():
        g,l,s=m.physical(u,w,off)
        return dict(error2=float(m.loss(u,w).sqrt()),gap=float(g),Lambda=float(l),lambda_eff=float(s.sqrt()) if s>0 else None,weight_norm=float(w.norm()))
