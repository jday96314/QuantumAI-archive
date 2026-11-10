"""Continuous/decoded PSD rank scans, retaining exact checkpoints and bounds."""
import argparse
from core import *
from solve import solve
from angle_gauge_probe import encode

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--tag',required=True)
    p.add_argument('--quantized',action='store_true');p.add_argument('--ranks',type=int,nargs='+',default=[88,128,144,160,175,192,224,256])
    p.add_argument('--steps',type=int,default=600);p.add_argument('--full-steps',type=int,default=2400)
    p.add_argument('--warm');a=p.parse_args()
    z=dict(np.load(a.input));u=z['U'];tree=z.get('tree_edges',np.load(INC)['tree_edges']);words=None
    if a.quantized:
        words=encode(u.reshape(-1,76),tree,8,True);u=reconstruct_u_vectors(words/128,tree).reshape(u.shape)
    m=Covariance(u);S=m.from_weights(z['w']) if 'w' in z else None
    if 'S' in z and not a.quantized:S=tensor(z['S'])
    S,hf=solve(m,S,a.full_steps,None,'fista',200)
    full=m.diagnostics(S,True);kwargs={} if words is None else dict(angle_words=words,brot=8)
    np.savez(HERE/f'parameters/{a.tag}_psd.npz',U=u,w=m.weights(S).cpu().numpy(),S=S.cpu().numpy(),tree_edges=tree,**kwargs)
    records=[]
    report(HERE/f'results/{a.tag}_capacity.json',dict(input=a.input,quantized=a.quantized,full_PSD=full,full_history=hf,ranks=records))
    warmz=dict(np.load(a.warm)) if a.warm else z
    warmS=m.from_weights(warmz['w']) if 'w' in warmz else None
    for C in a.ranks:
        if C>m.B:continue
        St=project(S,C)
        if float((St-S).norm())<1e-10:Sc,hc=St,[]
        else:Sc,hc=solve(m,S,a.steps,C,'coordinate',200)
        # Compare a continuation of the existing rank-C basin. Truncating the
        # convex optimum can be a worse nonconvex initialization.
        if warmS is not None:
            Sw,hw=solve(m,warmS,a.steps,C,'coordinate',200)
            if float(m.objective(Sw))<float(m.objective(Sc)):Sc,hc=Sw,hw
        diag=m.diagnostics(Sc,True,C);records.append(dict(C=C,**diag))
        np.savez(HERE/f'parameters/{a.tag}_c{C}.npz',U=u,w=m.weights(Sc,C).cpu().numpy(),S=Sc.cpu().numpy(),tree_edges=tree,**kwargs)
        report(HERE/f'results/{a.tag}_c{C}.json',dict(rank=C,diagnostics=diag,history=hc))
        report(HERE/f'results/{a.tag}_capacity.json',dict(input=a.input,quantized=a.quantized,full_PSD=full,full_history=hf,ranks=records))
if __name__=='__main__':main()
