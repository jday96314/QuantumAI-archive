"""Inputs and symmetric coordinates, independent of mutable sibling scripts."""
from pathlib import Path
import json
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
SOURCE=ROOT/'initial_experiments/4/codex/candidates/norm_cycle2_signed/submission.npz'
SPECTRUM=ROOT/'initial_experiments/3/parameters/reference_spectrum.npz'
REFERENCE=ROOT/'initial_experiments/3/parameters/chem_reference.npz'

def load():
    from dfthcsos.scoreinfra.dfthc_math import reconstruct_u_vectors
    z=dict(np.load(SOURCE));sp=dict(np.load(SPECTRUM))
    z['U']=reconstruct_u_vectors(z['angle_words'][76:]/128,z['tree_edges']).reshape(1,719,76)
    z['w']=z['W'].copy();z['offset']=z['H2_I']-z['w'].sum(1)
    return z,sp

def feature(u,sp):return u[...,sp['ii']]*u[...,sp['jj']]*sp['fac']

def unpack(f,sp):
    out=np.zeros((*f.shape[:-1],76,76))
    out[...,sp['ii'],sp['jj']]=f/sp['fac']
    out[...,sp['jj'],sp['ii']]=f/sp['fac']
    return out

def gram(sp):return (sp['eigenvectors']*sp['eigenvalues'])@sp['eigenvectors'].T

def report(name,data):
    p=HERE/'results'/name;p.parent.mkdir(exist_ok=True)
    p.write_text(json.dumps(data,indent=2)+'\n')

def emit(data):print(json.dumps(data),flush=True)
