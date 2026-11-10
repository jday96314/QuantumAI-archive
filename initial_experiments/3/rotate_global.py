"""Apply sampled Jacobi column rotations, preserving the represented tensor."""
import argparse,subprocess
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',required=True)
p.add_argument('--steps',type=int,default=500000)
p.add_argument('--method',choices=['sample','exact'],default='sample');a=p.parse_args();HERE=Path(__file__).resolve().parent
z=dict(np.load(a.fit));w=z['w'];W=z['W'];R,B,C=w.shape
rotw=[];rotW=[]
for r in range(R):
    A=np.concatenate([w[r],(W[r]+w[r].sum(0))[None,:]],axis=0).T.copy()
    src=Path(a.out).with_suffix(f'.{r}.in.bin');dst=Path(a.out).with_suffix(f'.{r}.out.bin');A.tofile(src)
    executable='rotate_exact' if a.method=='exact' else 'rotate_columns'
    subprocess.run([str(HERE/executable),str(src),str(dst),str(C),str(B+1),str(a.steps)],check=True)
    out=np.fromfile(dst,dtype=float).reshape(C,B+1).T
    assert np.allclose(A.T@A,out@out.T,rtol=1e-8,atol=1e-9)
    rotw.append(out[:B]);rotW.append(out[B]-out[:B].sum(0))
z['w']=np.array(rotw);z['W']=np.array(rotW);np.savez(a.out,**z)
