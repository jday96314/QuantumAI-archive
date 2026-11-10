"""Test exact gradient of captured Frobenius energy with respect to discrete angles theta."""

import time, torch, numpy as np

device = torch.device("cuda")
sp = np.load("../../3/parameters/reference_spectrum.npz")
ii = torch.tensor(sp["ii"], device=device)
jj = torch.tensor(sp["jj"], device=device)
fac = torch.tensor(sp["fac"], dtype=torch.float64, device=device)
ev = torch.tensor(sp["eigenvalues"], dtype=torch.float64, device=device)
V = torch.tensor(sp["eigenvectors"], dtype=torch.float64, device=device)
G = (V * torch.clamp(ev, min=0.0)) @ V.T
total_norm_sq = torch.sum(torch.clamp(ev, min=0.0)**2).item()
tree = np.load("../../3/candidates/chem35c45_rot.npz")["tree_edges"]

stages = []
for j, (pa, ch) in enumerate(tree):
    pa, ch = int(pa), int(ch)
    if not stages or any(pa in (x, y) or ch in (x, y) for _, x, y in stages[-1]):
        stages.append([])
    stages[-1].append((j, pa, ch))
stages = [tuple(torch.tensor(x, device=device) for x in zip(*s)) for s in stages]

B = 400
bits = 8
alpha = np.pi / (1 << (bits - 1))

# Initialize theta from grow650 top 400 atoms
z = np.load("../../3/candidates/grow650_q8_c350.npz")
words = z["angle_words"][:B]
theta = torch.tensor(words, dtype=torch.float64, device=device, requires_grad=True)

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

un = decode(theta)
P = un[:, ii] * un[:, jj] * fac
S = P @ P.T + 1e-6 * torch.eye(B, dtype=torch.float64, device=device)
L = torch.linalg.cholesky(S)
H = P @ G @ P.T
Y = torch.linalg.solve_triangular(L, H, upper=False)
Y = torch.linalg.solve_triangular(L, Y.T, upper=False).T
inside_sq = torch.sum(Y**2)
loss = -inside_sq # maximize inside_sq (minimize floor!)

loss.backward()

print("Gradient computation successful!")
print(f"Loss = {loss.item():.4f}, current floor = {np.sqrt(max(0.0, total_norm_sq - inside_sq.item()))*1000:.3f} mHa")
print(f"Grad norm on theta: {torch.norm(theta.grad).item():.6f}")
