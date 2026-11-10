"""Read-only round-4 selected inputs and physical symmetric feature coordinates."""
from pathlib import Path
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ROUND4 = ROOT / 'initial_experiments/4/codex'
SOURCE = ROUND4 / 'candidates/norm_cycle2_signed/submission.npz'
SPECTRUM = ROOT / 'initial_experiments/3/parameters/reference_spectrum.npz'
REFERENCE = ROOT / 'initial_experiments/3/parameters/chem_reference.npz'

def load():
    from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
    z = dict(np.load(SOURCE))
    z['U'] = reconstruct_u_vectors(z['angle_words'][76:] / 128, z['tree_edges']).reshape(1,719,76)
    z['w'] = z['W'].copy()
    z['offset'] = z['H2_I'] - z['w'].sum(1)
    return z, dict(np.load(SPECTRUM))

def feature(u, sp):
    return u[...,sp['ii']] * u[...,sp['jj']] * sp['fac']

def unpack(f, sp):
    out=np.zeros((*f.shape[:-1],76,76))
    out[...,sp['ii'],sp['jj']]=f/sp['fac']
    out[...,sp['jj'],sp['ii']]=f/sp['fac']
    return out

def report(path, data):
    Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(json.dumps(data,allow_nan=False),flush=True)
