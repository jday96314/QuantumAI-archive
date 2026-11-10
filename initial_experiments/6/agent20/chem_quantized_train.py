"""Train exact integer angle words from the accepted physical dictionary.

The forward pass decodes rounded integer words. Straight-through derivatives
propose joint word/weight/offset moves; float64 native-formula checkpoints
decide which research states to retain. No surrogate score is accepted.
"""
import argparse
import copy
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import torch
from audit_core import Metric
from chem_rank_screen import HERE,BANK,quantities


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--steps',type=int,default=2000)
    parser.add_argument('--lambda-cap',type=float,default=49.)
    parser.add_argument('--tag',default='quantized_control')
    parser.add_argument('--angle-lr',type=float,default=.015)
    parser.add_argument('--weight-lr',type=float,default=3e-5)
    parser.add_argument('--offset-lr',type=float,default=3e-5)
    parser.add_argument('--warmup',type=int,default=0)
    parser.add_argument('--project-gap',action='store_true')
    args=parser.parse_args()
    (HERE/'parameters'/f'{args.tag}_run.json').write_text(json.dumps(dict(
        arguments=vars(args),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        baseline_sha256=hashlib.sha256(BANK.read_bytes()).hexdigest(),
        precision='float32 training; exact decoded float64 physical checkpoints'),indent=2)+'\n')
    torch.set_num_threads(1);torch.manual_seed(20261025)
    z=dict(np.load(BANK));exact=Metric(z['tree_edges'])
    m=copy.copy(exact);m.d=dict(device='cuda',dtype=torch.float32)
    for name in ('fac','G','h0'):setattr(m,name,getattr(exact,name).float())
    original=z['angle_words'][76:].reshape(1,719,75)
    angles=torch.nn.Parameter(m.tensor(original))
    weights=torch.nn.Parameter(m.tensor(z['W']))
    offsets=torch.nn.Parameter(m.tensor(z['H2_I']-z['W'].sum(1)))
    optimizer=torch.optim.Adam([dict(params=[angles],lr=args.angle_lr),dict(params=[weights],lr=args.weight_lr),
                                dict(params=[offsets],lr=args.offset_lr)])
    started=time.time();history=[];best=float('inf')
    for step in range(args.steps+1):
        if step%20==0 or step==args.steps:
            with torch.no_grad():
                words=angles.round().cpu().numpy().astype(np.int64)%256
                u=exact.decode(exact.tensor(words));w=exact.tensor(weights.detach().cpu().numpy())
                off=exact.tensor(offsets.detach().cpu().numpy())
                if args.project_gap and step:
                    exact.project_gap(u,w,off,target=3.91515)
                    offsets.copy_(off.float())
                row=dict(step=step,seconds=time.time()-started,
                         changed_words=int(np.count_nonzero(words!=original)),**quantities(exact,u,w,off))
                history.append(row);print(json.dumps(row),flush=True)
                if row['lambda_eff']<=args.lambda_cap and row['gap']>=3.91504 and row['error2']<best:
                    best=row['error2']
                    np.savez(HERE/'parameters'/f'{args.tag}_best.npz',U=u.cpu().numpy(),w=w.cpu().numpy(),
                             W=off.cpu().numpy(),tree_edges=z['tree_edges'],angle_words=words.reshape(-1,75),brot=8)
                    (HERE/'results'/f'{args.tag}_best.json').write_text(json.dumps(row,indent=2)+'\n')
                (HERE/'results'/f'{args.tag}_history.json').write_text(json.dumps(history,indent=2)+'\n')
            if step==args.steps:break
        optimizer.zero_grad(set_to_none=True)
        quantized=angles+(angles.round()-angles).detach()
        u=m.decode(quantized)
        error2=m.loss(u,weights)
        gap,lam,effective=m.physical(u,weights,offsets)
        lambda_target=(args.lambda_cap**2+3.91515**2)/(2*3.91515)
        normalization=torch.relu(lam/lambda_target-1) if args.project_gap else torch.relu(effective/(args.lambda_cap**2)-1)
        gap_shortfall=torch.relu((3.9151-gap)/3.9151)
        objective=error2/.0004+1000*normalization.square()
        if not args.project_gap:objective+=100000*gap_shortfall.square()
        objective.backward()
        if step<args.warmup:angles.grad.zero_()
        torch.nn.utils.clip_grad_norm_([angles,weights,offsets],100.)
        optimizer.step()
        with torch.no_grad():
            angles.clamp_(m.tensor(original)-2,m.tensor(original)+2)
