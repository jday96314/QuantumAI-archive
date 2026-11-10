"""Necessary tensor-error floor for one exact-word dictionary, any weights.

Every allowed factor lies in the span of its rank-one projector features.
Orthogonal projection of the target tensor to that union span bounds all
block-PSD and rank-constrained fits from below. It is not a shape-wide bound.
"""
import argparse,json,time,hashlib
from pathlib import Path
from audit_core import *
HERE=Path(__file__).resolve().parent
@torch.no_grad()
def main():
 p=argparse.ArgumentParser();p.add_argument('--init',required=True);p.add_argument('--tag',required=True);a=p.parse_args();torch.set_num_threads(1);start=time.time();z,u,w,tree,words=load(a.init);m=Metric(tree)
 if words is None:words=encode(u.reshape(-1,76),tree,8,True)
 uq=m.decode(m.tensor(words).reshape(w.shape[0],w.shape[1],75));P=m.feats(uq).flatten(0,1);_,sv,V=torch.linalg.svd(P,full_matrices=False);rank=int((sv>sv[0]*1e-10).sum());Q=V[:rank];small=Q@m.G@Q.T;remainder=(m.G.square().sum()-small.square().sum()).clamp(min=0).sqrt()
 r=dict(status='FIXED_WORD_UNION_SPAN_NECESSARY_BOUND',init=a.init,shape=list(w.shape),dictionary_rank=rank,smallest_singular_value=float(sv[-1]),union_span_error_floor=float(remainder),native_current_error=float(m.loss(uq,m.tensor(w)).sqrt()),seconds=time.time()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),parameter_sha256=hashlib.sha256(Path(a.init).read_bytes()).hexdigest(),scope='orthogonal projection bound for this exact rounded dictionary with arbitrary signed full-rank union weights; no exclusion of another geometry or whole shape; no native circuit')
 (HERE/'results'/f'{a.tag}.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':main()
