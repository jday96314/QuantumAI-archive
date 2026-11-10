"""Refine weighted ensemble truncations near the accuracy boundary."""
import numpy as np
z=dict(np.load('initial_experiments/candidates/chemistry-weighted-g0.1-c32.npz'))
for C in [29,30]:
 z=dict(np.load('initial_experiments/candidates/chemistry-weighted-g0.1-c32.npz'))
 zz={**z,'W':z['W'][:,:,:C],'H2_I':z['H2_I'][:,:C]}
 start=f'chemistry-c{C}-weighted-start.npz';out=f'chemistry-c{C}-weighted.npz'
 np.savez('initial_experiments/candidates/'+start,**zz)
 source=open('initial_experiments/chemistry_refine_rotations.py').read().replace('chemistry-rotated.npz',start).replace('chemistry-refined.npz',out)
 source=source.replace('range(100)','range(100)').replace('scale=[.10,.04,.015,.006,.002,.0008,.0003,.0001][min(sweep//12,7)]','scale=[.7854,.3,.10,.04,.015,.006,.002,.0008,.0003,.0001][min(sweep//10,9)]')
 exec(compile(source,f'chemistry_c{C}_refine','exec'),{'__name__':f'chemistry_c{C}_refine'})
 print('COMPLETE C',C,flush=True)
