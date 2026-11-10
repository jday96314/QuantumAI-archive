"""Deterministic-seeded annealing for simultaneous integer angle proposals."""
import ctypes
from pathlib import Path
import numpy as np
lib=ctypes.CDLL(str(Path(__file__).resolve().parent/'qubo_anneal.so'))
array=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
lib.solve_qubo.argtypes=[ctypes.c_int,array,array,ctypes.c_int,ctypes.c_int,ctypes.c_uint64,np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')]
lib.solve_qubo.restype=ctypes.c_double
def solve(H,b,restarts=12,cap=64,seed=71842):
 H=np.ascontiguousarray(.5*(H+H.T),dtype=np.float64);np.fill_diagonal(H,0)
 b=np.ascontiguousarray(b,dtype=np.float64);x=np.zeros(len(b),dtype=np.int32)
 value=lib.solve_qubo(len(b),H,b,restarts,cap,seed,x)
 # Independently recompute the proposed energy; do not trust a drifting field.
 exact=float(b@x+.5*x@H@x)
 scale=max(float(np.max(np.abs(H))),float(np.max(np.abs(b))),1e-30)
 atol=max(1e-15,64*np.finfo(float).eps*len(b)**2*scale)
 assert np.isclose(value,exact,atol=atol,rtol=1e-5),(value,exact)
 return np.flatnonzero(x),exact
