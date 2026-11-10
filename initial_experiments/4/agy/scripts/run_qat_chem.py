"""Quantization-Aware Training (QAT) at brot=8 on GPU.

Jointly optimizes discrete angle words and factorization weights w to cross err2 <= 19.95 mHa.
"""

import argparse, json, time, os, sys
from pathlib import Path

ROOT = Path("/media/shared/ExpansionDrive2/QuantumAI")
sys.path[:0] = [
    str(ROOT / "from_host/kickmix"),
    str(ROOT / "from_host/verifiers"),
    str(ROOT / "initial_experiments/3"),
]

import numpy as np
import torch
from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fit", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--C", type=int, default=350)
    p.add_argument("--brot", type=int, default=8)
    p.add_argument("--steps", type=int, default=2000)
    p.add_argument("--angle-lr", type=float, default=0.005)
    p.add_argument("--weight-lr", type=float, default=5e-5)
    a = p.parse_args()
    
    device = torch.device("cuda")
    z = dict(np.load(a.fit))
    tree = z["tree_edges"]
    brot = a.brot
    alpha = np.pi / (1 << (brot - 1))
    
    # Orbitals
    N = 76
    rows, cols = np.tril_indices(N)
    weights_tril = np.where(rows == cols, 1.0, np.sqrt(2.0))
    weights_t = torch.tensor(weights_tril, dtype=torch.float64, device=device)
    rows_t = torch.tensor(rows, dtype=torch.long, device=device)
    cols_t = torch.tensor(cols, dtype=torch.long, device=device)
    
    # Reference G
    sp = np.load(ROOT / "initial_experiments/3/parameters/reference_spectrum.npz")
    ev, V = sp["eigenvalues"], sp["eigenvectors"]
    G = (V * np.maximum(ev, 0)) @ V.T # (2926, 2926)
    G_t = torch.tensor(G, dtype=torch.float64, device=device)
    total_norm_sq = float(np.sum(np.maximum(ev, 0)**2))
    
    # Load or extract words and w
    words = z["angle_words"].reshape(-1, 75)
    B = len(words)
    R = 1
    C = a.C
    
    if "w" in z:
        w_init = z["w"].reshape(R, B, -1)[:, :, :C]
    elif "W" in z and z["W"].ndim == 3:
        w_init = z["W"][:, :, :C]
    elif "K" in z:
        K = z["K"]
        e_k, v_k = np.linalg.eigh(K)
        top_v = v_k[:, -C:]
        top_e = np.maximum(e_k[-C:], 0.0)
        w_init = (top_v * np.sqrt(top_e))[None, :, :]
    else:
        raise ValueError("Cannot find w, W, or K in fit file.")
        
    print(f"Loaded candidate: R={R}, B={B}, C={C}, brot={brot}")
    
    # Build Givens stages
    stages = []
    for t, (pa, ch) in enumerate(tree):
        pa, ch = int(pa), int(ch)
        if not stages or any(pa in (p0, c0) or ch in (p0, c0) for _, p0, c0 in stages[-1]):
            stages.append([])
        stages[-1].append((t, pa, ch))
    stages = [tuple(torch.tensor(v, device=device) for v in zip(*st)) for st in stages]
    
    theta = torch.tensor(words, device=device, dtype=torch.float64, requires_grad=True)
    w = torch.tensor(w_init, device=device, dtype=torch.float64, requires_grad=True)
    
    eye0 = torch.zeros(R * B, N, device=device, dtype=torch.float64)
    eye0[:, 0] = 1.0
    
    def decode(t):
        rounded = t + (t.round() - t).detach()
        angles = rounded * alpha
        u = eye0.clone()
        for ts, ps, cs in stages:
            ang = angles[:, ts]
            co, sn = ang.cos(), ang.sin()
            up, uc = u[:, ps], u[:, cs]
            nxt = u.clone()
            nxt[:, ps] = co * up - sn * uc
            nxt[:, cs] = sn * up + co * uc
            u = nxt
        return u
        
    def loss():
        un = decode(theta)
        P = (un[:, rows_t] * un[:, cols_t]) * weights_t[None, :] # (R*B, 2926)
        P_r = P.reshape(R, B, 2926)
        X = torch.bmm(w.transpose(1, 2), P_r).reshape(R * C, 2926)
        diff = X.T @ X - G_t
        return torch.sum(diff * diff)
        
    # Evaluate initial discrete floor
    with torch.no_grad():
        e2_0 = loss()
        err2_0 = float(e2_0.sqrt())
        print(f"Initial exact discrete err2: {err2_0*1000:.4f} mHa")
        
    opt = torch.optim.Adam([
        {"params": [theta], "lr": a.angle_lr},
        {"params": [w], "lr": a.weight_lr},
    ])
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=a.steps, eta_min=1e-6)
    
    t0 = time.time()
    best = err2_0
    
    for s in range(a.steps):
        opt.zero_grad()
        e2 = loss()
        val = float(e2.detach().sqrt())
        
        if val < best:
            best = val
            # Snapshot best solution
            zz = dict(z)
            words_round = theta.detach().round().cpu().numpy().astype(np.int64) % (1 << brot)
            u_eval = decode(theta).detach().cpu().numpy().reshape(R, B, N)
            w_eval = w.detach().cpu().numpy()
            zz.update(
                U=u_eval,
                w=w_eval,
                angle_words=words_round,
                brot=brot,
                err2=val,
            )
            np.savez(a.out, **zz)
            
        if val * 1000 <= 19.95:
            print(f"\n*** SUCCESS: Reached discrete err2 = {val*1000:.4f} mHa < 19.95 mHa at step {s}! ***")
            break
            
        e2.backward()
        opt.step()
        sched.step()
        
        if (s + 1) % 50 == 0:
            print(f"Step {s+1:4d} | discrete err2 = {val*1000:6.3f} mHa | best = {best*1000:6.3f} mHa | Elapsed: {time.time()-t0:.1f}s")
            
    print(f"\nQAT complete in {time.time()-t0:.1f}s. Best discrete err2 = {best*1000:.4f} mHa.")
    print(f"Saved snapshot to {a.out}.")

if __name__ == "__main__":
    main()
