"""Measure sign mixing and exact factor metrics without changing a candidate."""
import argparse,hashlib,io,json
from pathlib import Path
import numpy as np

def measure(path):
 root=Path(__file__).resolve().parent
 with np.load(io.BytesIO(Path(path).read_bytes())) as z:z={k:z[k] for k in z.files}
 u,w,W=z['U'],z['w'],z['W'];sp=np.load(root/'parameters/reference_spectrum.npz');ref=np.load(root/'parameters/chem_reference.npz')
 ii,jj,fac=sp['ii'],sp['jj'],sp['fac'];P=u[:,:,ii]*u[:,:,jj]*fac
 X=np.einsum('rbc,rbp->rcp',w,P).reshape(-1,len(ii))
 G=(sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T
 er=float(np.linalg.norm(X.T@X-G));l1=float((abs(w).sum(1)**2).sum());tr=float((w.sum(1)**2).sum())
 packed=np.einsum('rc,rcp->p',W,X.reshape(w.shape[0],w.shape[2],-1))/fac
 correction=np.zeros((76,76));correction[ii,jj]=packed;correction+=correction.T-np.diag(correction.diagonal())
 h=ref['h1_thc']-.5*int(ref['electrons'])*ref['S']+float(ref['beta_1'])*np.eye(76)-correction
 e=np.linalg.eigvalsh(h);lh1=float(abs(e).sum())
 constant=-22140.27335-float(ref['ecore'])+int(ref['electrons'])*float(ref['beta_1'])
 gap=constant+lh1-e.sum()+.5*(W*W).sum()
 h2=float(.25*((abs(w).sum(1)+abs(W+w.sum(1)))**2).sum());lam=lh1+h2
 return dict(path=str(path),shape=list(w.shape),brot=int(z.get('brot',8)),mixing=l1-tr,l1_squared=l1,trace_squared=tr,error_2=er,gap=float(gap),Lambda=lam,lambda_eff_model=float(np.sqrt(gap*(2*lam-gap))),h1_norm=lh1,h2_norm=h2,U_sha256=hashlib.sha256(u.tobytes()).hexdigest(),tree_sha256=hashlib.sha256(z['tree_edges'].tobytes()).hexdigest(),angle_words_sha256=hashlib.sha256(z['angle_words'].tobytes()).hexdigest())

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('fits',nargs='+');p.add_argument('--output',type=Path);a=p.parse_args()
 out=[measure(f) for f in a.fits];text=json.dumps(out,indent=2)+'\n'
 if a.output:a.output.write_text(text)
 print(text,end='')
