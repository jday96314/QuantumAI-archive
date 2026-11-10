import numpy as np,sys
from dfthcsos.verify import load_inputs, evaluate_words
inp=load_inputs('parameters/chemistry.npz')
d=np.load('parameters/chemistry.npz');w=d['angle_words'].astype(np.int64);tree=d['tree_edges']
for brot in [12,11,10,9]:
    sh=12-brot
    ww=((w+(1<<sh>>1))>>sh)%(1<<brot) if sh else w
    m=evaluate_words(inp,ww,tree,brot)
    print(brot,{k:(round(v,6) if isinstance(v,float) else v) for k,v in m.items()},flush=True)
