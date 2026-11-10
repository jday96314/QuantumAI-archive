"""Test speed and gradient of differentiable Givens tree decoder on GPU."""

import time, torch, numpy as np

device = torch.device("cuda")
tree = np.load("../../3/candidates/chem35c45_rot.npz")["tree_edges"]

stages = []
for j, (pa, ch) in enumerate(tree):
    pa, ch = int(pa), int(ch)
    if not stages or any(pa in (x, y) or ch in (x, y) for _, x, y in stages[-1]):
        stages.append([])
    stages[-1].append((j, pa, ch))
stages = [tuple(torch.tensor(x, device=device) for x in zip(*s)) for s in stages]
print(f"Tree grouped into {len(stages)} parallel stages.")

B = 400
bits = 8
alpha = np.pi / (1 << (bits - 1))
theta = torch.randn(B, 75, dtype=torch.float64, device=device, requires_grad=True)

def decode(th):
    angles = (th + (th.round() - th).detach()) * alpha
    un = torch.zeros((len(th), 76), dtype=torch.float64, device=device)
    un[:, 0] = 1.0
    for ts, ps, cs in stages:
        ang = angles[:, ts]
        c, s = ang.cos(), ang.sin()
        up, uc = un[:, ps], un[:, cs]
        new = un.clone()
        new[:, ps] = c * up - s * uc
        new[:, cs] = s * up + c * uc
        un = new
    return un

# Benchmark forward and backward
torch.cuda.synchronize()
t0 = time.time()
for _ in range(100):
    un = decode(theta)
    loss = torch.sum(un**2)
    loss.backward()
    theta.grad.zero_()
torch.cuda.synchronize()
t1 = time.time()
print(f"100 forward+backward passes took {t1 - t0:.3f}s ({(t1 - t0)*10:.2f} ms per step)!")
