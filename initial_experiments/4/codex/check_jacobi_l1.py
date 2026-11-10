"""Compare the exact rotation search with a dense independent angle grid."""
import json,subprocess
from pathlib import Path
import numpy as np
from jacobi_l1 import encode,decode,root

rng=np.random.default_rng(719281)
theta=np.linspace(0,np.pi/2,32769);co=np.cos(theta);si=np.sin(theta)
max_gram_error=0.;cases=0
for n in [2,3,5,12,32]:
    for trial in range(12):
        A=rng.normal(size=(n,2))
        if trial%3==0:A[rng.random(A.shape)<.35]=0
        ip=root/'parameters/jacobi_check_input.bin';op=root/'parameters/jacobi_check_output.bin'
        encode(A,ip)
        subprocess.run([str(root/'jacobi_l1'),str(ip),str(op),'1','1','11271'],check=True,stdout=subprocess.DEVNULL)
        out=decode(op)
        xx=A[:,0,None]*co+A[:,1,None]*si
        yy=-A[:,0,None]*si+A[:,1,None]*co
        grid=.25*(abs(xx).sum(0)**2+abs(yy).sum(0)**2)
        got=float(.25*np.square(abs(out).sum(0)).sum())
        assert got<=float(grid.min())+1e-10,(n,trial,got,grid.min())
        err=float(abs(A@A.T-out@out.T).max());max_gram_error=max(max_gram_error,err)
        assert err<1e-12,(n,trial,err)
        cases+=1
result=dict(status='PASS',cases=cases,angles_per_brute_grid=len(theta),max_gram_error=max_gram_error)
(root/'results/jacobi_grid_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)
