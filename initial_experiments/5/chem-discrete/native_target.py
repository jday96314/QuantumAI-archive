"""Build the exact symmetric-pair target from unchanged native reference inputs."""
from pathlib import Path
import os,json
import numpy as np
from dfthcsos.verify import load_inputs
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def build():
    inc=ROOT/'initial_experiments/4/codex/candidates/norm_cycle2_signed/submission.npz'
    inp=load_inputs(inc);sp=np.load(ROOT/'initial_experiments/3/parameters/reference_spectrum.npz')
    ii,jj,fac=sp['ii'],sp['jj'],sp['fac'];assert np.allclose(fac,np.where(ii==jj,1,np.sqrt(2)),rtol=0,atol=1e-15)
    G=inp.ref['tei'][ii[:,None],jj[:,None],ii[None,:],jj[None,:]]*fac[:,None]*fac[None,:]
    identity=(ii==jj).astype(float);shift=inp.S[ii,jj]*fac
    G+=.5*(shift[:,None]*identity[None,:]+identity[:,None]*shift[None,:])
    path=HERE/'candidates/native_pair_target.npz';tmp=path.with_name(path.name+f'.{os.getpid()}.tmp')
    with tmp.open('wb') as f:np.savez(f,G=G,ii=ii,jj=jj,fac=fac)
    tmp.replace(path)
    h0=inp.ref['h1_thc']-.5*inp.ref['electrons']*inp.S+inp.beta_1*np.eye(76)
    constant=inp.energy_upper-inp.ref['ecore']+inp.ref['electrons']*inp.beta_1
    path=HERE/'candidates/native_constants.npz';tmp=path.with_name(path.name+f'.{os.getpid()}.tmp')
    with tmp.open('wb') as f:np.savez(f,h0=h0,constant=constant)
    tmp.replace(path)
    return G
if __name__=='__main__':
    G=build();print(json.dumps(dict(shape=list(G.shape),source='native checksum-validated float32 TEI and fixed symmetry shifts')))
