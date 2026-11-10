from pathlib import Path
import struct,numpy as np
from quadratic_schur_model import F,K,node_maps
from dqi.verify_submission_dqi import optimal_dqi_weights,default_num_accepted_values
nodes,kap,P,maps=node_maps(35);prob=optimal_dqi_weights(4095,35,12,default_num_accepted_values(12))**2
with Path('parameters/schur_census.bin').open('wb') as f:
 f.write(struct.pack('<II',70,K));f.write(np.asarray(prob,dtype='<f8').tobytes());f.write(np.asarray(F.ex[:4095],dtype='<u2').tobytes());f.write(np.asarray(F.lg,dtype='<u2').tobytes());f.write(np.asarray(kap,dtype='<u2').tobytes());f.write(np.asarray(maps,dtype='<u4').tobytes());f.write(np.asarray(P,dtype='<u2').tobytes())
