"""Exact plane-rotation search on the native checkpoint's column gauge.

The combined weight/identity coefficient Gram matrix is preserved. Existing
angle circuits can be reused unchanged if the full chemistry verifier passes.
"""
import argparse,json,struct,subprocess,time
from pathlib import Path
import numpy as np

root=Path(__file__).resolve().parent


def encode(A,path):
    Path(path).write_bytes(struct.pack('<II',*A.shape)+np.asarray(A,dtype='<f8').tobytes())


def decode(path):
    data=Path(path).read_bytes();n,c=struct.unpack('<II',data[:8])
    return np.frombuffer(data[8:],dtype='<f8').reshape(n,c).copy()


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',default='candidates/selectswap7_anf/submission.npz')
    p.add_argument('--output',default='parameters/jacobi_native.npz');p.add_argument('--passes',type=int,default=12)
    p.add_argument('--pairs',type=int,default=12000);p.add_argument('--seed',type=int,default=791291)
    a=p.parse_args();z=dict(np.load(a.input));w=z['W'];assert w.shape[0]==1
    A=np.concatenate([w[0],z['H2_I'][0,None,:]],axis=0)
    inputfile=root/'parameters/jacobi_input.bin';outputfile=root/'parameters/jacobi_output.bin'
    encode(A,inputfile);began=time.time()
    subprocess.run([str(root/'jacobi_l1'),str(inputfile),str(outputfile),str(a.passes),str(a.pairs),str(a.seed)],check=True)
    rot=decode(outputfile)
    gram=A@A.T;delta=rot@rot.T-gram
    check=dict(status='PASS' if np.max(abs(delta))<1e-9 else 'FAIL',gram_max_error=float(abs(delta).max()),
               lambda_h2_before=float(.25*np.square(abs(A).sum(0)).sum()),
               lambda_h2_after=float(.25*np.square(abs(rot).sum(0)).sum()),seconds=time.time()-began)
    assert check['status']=='PASS',check
    z['W']=rot[:-1][None,:,:];z['H2_I']=rot[-1][None,:]
    np.savez(root/a.output,**z)
    (root/'results/jacobi_gram_check.json').write_text(json.dumps(check,indent=2)+'\n')
    print(json.dumps(check),flush=True)


if __name__=='__main__':main()
