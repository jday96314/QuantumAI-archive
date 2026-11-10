"""Convex upper-gap offset probe and paths back to the verified lower branch.

For fixed physical factors, 2*Lambda-Egap is convex and piecewise linear in
offsets except for the convex positive-eigenvalue sum. Its quadratic terms
cancel exactly. Optimize it directly; evaluate the true product and native
floors independently. Research data only.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import torch


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fit', default='parameters/norm_lattice719_final_fit.npz')
    p.add_argument('--steps', type=int, default=10000)
    p.add_argument('--lr', type=float, default=.002)
    p.add_argument('--output', default='candidates/upper_gap_offsets.npz')
    a = p.parse_args()
    torch.set_num_threads(2)
    z = dict(np.load(a.fit))
    ref = np.load('../../3/parameters/chem_reference.npz')
    u, w, W = z['U'][0], z['w'][0], z['W'][0]
    d = dict(device='cuda', dtype=torch.float64)
    X = torch.tensor(np.einsum('bp,bq,bc->cpq', u, u, w, optimize=True), **d)
    off = torch.tensor(W, requires_grad=True, **d)
    trace = torch.tensor(w.sum(0), **d)
    l1 = torch.tensor(abs(w).sum(0), **d)
    h0 = torch.tensor(ref['h1_thc']-.5*int(ref['electrons'])*ref['S']+float(ref['beta_1'])*np.eye(76), **d)
    const = -22140.27335-float(ref['ecore'])+int(ref['electrons'])*float(ref['beta_1'])

    def values(v):
        h = h0 - torch.einsum('c,cpq->pq', v, X)
        eig = torch.linalg.eigvalsh(h)
        gap = const + 2 * torch.relu(-eig).sum() + .5*v.square().sum()
        upper = -const + 2 * torch.relu(eig).sum() + (.5*(trace.square()+l1.square())+v*trace+l1*(v+trace).abs()).sum()
        lam = .5*(gap+upper)
        return gap, upper, lam, gap*upper

    best = float('inf')
    best_upper = float('inf')
    upper_off = None
    initial = values(off)
    initial_score = float(initial[3].detach())
    opt = torch.optim.Adam([off], lr=a.lr)
    began = time.time()
    for step in range(a.steps):
        opt.zero_grad()
        vals = values(off)
        vals[1].backward()
        opt.step()
        with torch.no_grad():
            if step % 10:
                continue
            vals = values(off)
            g, upper, lam, prod = map(float, vals)
            if upper < best_upper:
                best_upper = upper
                upper_off = off.detach().clone()
            if g > 3.914898 and lam > 49 and upper > 0 and prod < min(best, initial_score):
                best = prod
                out = dict(z)
                out['W'] = off.detach().cpu().numpy()[None]
                for key in ['onebody_angle_words', 'onebody_brot', 'error_1_optimized']:
                    out.pop(key, None)
                np.savez(a.output, **out)
            if step % 500 == 0:
                print(json.dumps(dict(step=step, gap=g, upper_gap=upper, Lambda=lam,
                                      lambda_eff=np.sqrt(max(0,prod)), best_valid_product=best,
                                      seconds=time.time()-began)), flush=True)
    with torch.no_grad():
        rows = []
        original = torch.tensor(W, **d)
        for t in np.linspace(0, 1.5, 601):
            g, upper, lam, prod = map(float, values(original+t*(upper_off-original)))
            rows.append(dict(t=float(t), gap=g, upper_gap=upper, Lambda=lam, product=prod))
        valid = [r for r in rows if r['gap']>3.914898 and r['upper_gap']>0 and r['Lambda']>49]
        record = dict(status='FINISHED', scope='convex fixed-factor upper-gap probe, not native verification',
                      initial=dict(zip(['gap','upper_gap','Lambda','product'],map(lambda x:float(x.detach()),initial))),
                      best_upper_gap=best_upper, best_valid_product=best,
                      best_path_point=min(valid,key=lambda x:x['product']), rows=rows)
        Path('results/'+Path(a.output).stem+'_probe.json').write_text(json.dumps(record,indent=2)+'\n')
        print(json.dumps({k:v for k,v in record.items() if k!='rows'}),flush=True)


if __name__ == '__main__':
    main()
