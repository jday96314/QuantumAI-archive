"""Exercise exact integer-move acceptance against direct tensor reconstruction."""
import numpy as np
import torch
from common import HERE,SPECTRUM,report
from angle_gauge_probe import encode
from quantized_screen import integer_search

def main():
    torch.set_num_threads(2)
    z=dict(np.load(HERE/'parameters/r4b180c88_pca_adam.npz'))
    sp=np.load(SPECTRUM);d=dict(device='cuda',dtype=torch.float64)
    G=torch.tensor((sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T,**d)
    words=encode(z['U'].reshape(-1,76),z['tree_edges'],7,True)
    _,control=integer_search(words,torch.tensor(z['w'],**d),G,z['tree_edges'],7,
                             sp['ii'],sp['jj'],torch.tensor(sp['fac'],**d),4)
    assert control['accepted_moves']>0
    report(HERE/'results/integer_metric_control.json',dict(status='PASS',
        scope='Algebraic acceptance control on actual fit; cached target is sufficient for this identity test, not a physical candidate score.',
        assertion='Every accepted contracted tensor update agrees with independently rebuilt full tensor squared norm within 1e-8 and strictly improves it.',control=control))

if __name__=='__main__':main()
