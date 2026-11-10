"""Export a quantized global dictionary and build its ANGLE circuit family."""
import argparse,json
from pathlib import Path
import numpy as np
from compress_chem import encode
from chemistry_qrom import build
from dfthcsos.example_submission.angle_k_reference_circuit import batch_deltas
p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--out',type=Path,required=True)
p.add_argument('--batch',type=int);a=p.parse_args();HERE=Path(__file__).resolve().parent
z=dict(np.load(a.fit));ref=np.load(Path('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/parameters')/'chem_reference.npz');base=dict(np.load(Path('/media/shared/ExpansionDrive2/QuantumAI/initial_experiments/3/candidates')/'chem35c45_rot.npz'))
u,w,W=z['U'],z['w'],z['W'];R,B,C=w.shape;n=76;bits=int(z.get('brot',13))
h0=ref['h1_thc']-.5*int(ref['electrons'])*ref['S']+float(ref['beta_1'])*np.eye(n)
X=np.einsum('rbp,rbq,rbc->rcpq',u,u,w,optimize=True)
h=h0-np.einsum('rc,rcpq->pq',W,X,optimize=True)
tree=z.get('tree_edges',base['tree_edges'])
_,v=np.linalg.eigh(.5*(h+h.T));one=encode(v.T,tree,bits)
if 'onebody_angle_words' in z:
    assert int(z['onebody_brot'])==bits
    one=z['onebody_angle_words']
two=z.get('angle_words')
if two is None:two=encode(u.reshape(-1,n),tree,bits)
words=np.r_[one,np.asarray(two).reshape(R*B,75)]
out=dict(U=u,W=w,H2_I=W+w.sum(1),H1_I=base['H1_I'],Hsym_lower=base['Hsym_lower'],tree_edges=tree,angle_words=words)
a.out.mkdir(parents=True,exist_ok=True);np.savez(a.out/'submission.npz',**out)
address=(n+R*B-1).bit_length();batch=a.batch or max(1,(128-address-8)//bits)
for k,delta in enumerate(batch_deltas(words,batch)):
    c=build(delta,bits);assert c.num_qubits<=128,(k,c.num_qubits,batch,bits)
    c.to_file(str(a.out/f'angle_k_{k}.kmb'),format='kmb')
    print(json.dumps(dict(k=k,brot=bits,batch=batch,qubits=c.num_qubits,toffoli=c.operation_statistics().num_toffoli_like)),flush=True)
