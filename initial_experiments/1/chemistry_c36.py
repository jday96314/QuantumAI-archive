"""Compress ensemble rank, then optimize its orthogonal coordinates."""
import numpy as np
z=dict(np.load('from_host/sample_submission/dfthcsos/submission.npz'))
a=np.concatenate([z['W'],z['H2_I'][:,None,:]],axis=1).astype(float)
u,s,v=np.linalg.svd(a,full_matrices=False)
a=u[:,:,:36]*s[:,None,:36]
np.savez('initial_experiments/candidates/chemistry-c36-start.npz',**{**z,'W':a[:,:60,:],'H2_I':a[:,60,:]})
source=open('initial_experiments/chemistry_refine_rotations.py').read()
source=source.replace('chemistry-rotated.npz','chemistry-c36-start.npz').replace('chemistry-refined.npz','chemistry-c36-refined.npz').replace('range(100)','range(100)').replace("scale=[.10,.04,.015,.006,.002,.0008,.0003,.0001][min(sweep//12,7)]","scale=[.7854,.3,.10,.04,.015,.006,.002,.0008,.0003,.0001][min(sweep//10,9)]")
exec(compile(source,'chemistry_c36_refine','exec'))
