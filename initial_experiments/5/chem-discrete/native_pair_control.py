"""Compare the full cached pair target with native float32 reference semantics.

Re-evaluate the saved PSD dual certificate directly against the native target,
so its capacity conclusion does not depend on a cached eigendecomposition.
"""
import json,hashlib
from core import *
from dfthcsos.verify import load_inputs
inp=load_inputs(INC/'submission.npz');sp=np.load(SP);ii,jj,fac=sp['ii'],sp['jj'],sp['fac']
G=inp.ref['tei'][ii[:,None],jj[:,None],ii[None,:],jj[None,:]]*fac[:,None]*fac[None,:]
identity=(ii==jj).astype(float);shift=inp.S[ii,jj]*fac
G+=.5*(shift[:,None]*identity[None,:]+identity[:,None]*shift[None,:])
cached=(sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T;delta=G-cached
cert=np.load(HERE/'candidates/admm_fullPSD_B479_polish_dual.npz');D=cert['D'];words=cert['angle_words'];tree=cert['tree_edges'];u=reconstruct_u_vectors(words/128,tree).reshape(2,479,76)
P=u[:,:,ii]*u[:,:,jj]*fac;maxev=[]
for pp in P:
    Q=np.linalg.qr(pp.T)[0];maxev.append(float(np.linalg.eigvalsh(Q.T@D@Q)[-1]))
dual2=float(2*np.sum(G*D)-np.square(D).sum());row=dict(cached_target_Frobenius_difference=float(np.linalg.norm(delta)),cached_target_max_difference=float(abs(delta).max()),native_dual_lower_bound_squared=dual2,native_dual_lower_bound=float(np.sqrt(dual2)),dual_block_max_eigenvalues=maxev,excludes_native_20mHa=bool(dual2>.020**2 and max(maxev)<0),reference_spectrum_sha256=hashlib.sha256(SP.read_bytes()).hexdigest())
assert row['excludes_native_20mHa'],row
(HERE/'results/native_pair_control.json').write_text(json.dumps(row,indent=2));print(json.dumps(row))
