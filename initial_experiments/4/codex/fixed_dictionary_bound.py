"""Gauge-independent normalization certificate for a fixed projector dictionary.

For invertible feature matrix P and coefficient covariance H=P^-T G P^-1,
any sign matrix S with entries in [-1,1] certifies sum(abs(H))>=<S,H>.
The triangle inequality also gives sum_c ||w_c||_1^2>=sum(abs(H)).
Combining the native gap floor and physical Frobenius error gate gives a
numerical lower bound independent of factor-column rotations and rank.
"""
import json,time
from pathlib import Path
import numpy as np
import scipy.sparse as sparse
from dfthcsos.verify import heuristic_costs


def main():
 begin=time.time();z=np.load('../../3/parameters/reference_spectrum.npz');ii,jj=z['ii'],z['jj'];D=len(ii);n=76
 G=(z['eigenvectors']*z['eigenvalues'])@z['eigenvectors'].T;where={(int(i),int(j)):a for a,(i,j) in enumerate(zip(ii,jj))};rows=[];cols=[];vals=[]
 # B maps a physical packed symmetric factor to one-hot / pair coefficients.
 for p,(i,j) in enumerate(zip(ii,jj)):
  if i==j:rows.append(p);cols.append(p);vals.append(1.)
  else:
   rows.extend([p,where[(int(i),int(i))],where[(int(j),int(j))]]);cols.extend([p,p,p]);vals.extend([np.sqrt(2),-1/np.sqrt(2),-1/np.sqrt(2)])
 B=sparse.csr_matrix((vals,(rows,cols)),shape=(D,D));pr=[];pc=[];pv=[]
 for p,(i,j) in enumerate(zip(ii,jj)):
  if i==j:pr.append(p);pc.append(p);pv.append(1.)
  else:pr.extend([p,p,p]);pc.extend([p,where[(int(i),int(i))],where[(int(j),int(j))]]);pv.extend([1/np.sqrt(2),.5,.5])
 P=sparse.csr_matrix((pv,(pr,pc)),shape=(D,D));residual=P.T@B-sparse.eye(D);roundtrip=float(abs(residual.data).max()) if residual.nnz else 0.;assert roundtrip<1e-14,roundtrip
 H=B@G@B.T;S=np.sign(H);diag=(ii==jj).astype(float);A=(B.T@S)@B-np.outer(diag,diag)
 physical_epsilon=.020;linear=float((A*G).sum());penalty=physical_epsilon*float(np.linalg.norm(A));fixed_c=-559.1324731485038;u=-fixed_c+.5*(linear-penalty);gap_floor=3.914898;lam=np.sqrt(max(0,gap_floor*u))
 evidence=dict(status='NUMERICAL_CERTIFICATE',scope='all column gauges and ranks for the fixed full one-hot/pair dictionary, physical Frobenius error <=20mHa; floating-point evaluation of an analytic dual certificate',B=D,physical_epsilon=physical_epsilon,linear_reference=linear,dual_error_penalty=penalty,certified_u_lower_estimate=u,lambda_eff_lower_estimate=float(lam),native_gap_floor=gap_floor,formula='lambda_eff^2 > gap_floor * [-c + 0.5*(<A,Gref>-epsilon*||A||F)], A=B.T*sign(B*Gref*B.T)*B-d*d.T',resource_lower_estimates=[])
 for C in [276,300,350]:
  heur=heuristic_costs(n,1,D,C,3);evidence['resource_lower_estimates'].append(dict(C=C,heuristic_T=heur['toffoli'],heuristic_Q=heur['qubits'],cost_lower_without_angles=float(lam*heur['toffoli']*heur['qubits'])))
 evidence['feature_inverse_roundtrip_max_error']=roundtrip;evidence['seconds']=time.time()-begin;Path('results/full_pair_dictionary_bound.json').write_text(json.dumps(evidence,indent=2)+'\n');print(json.dumps(evidence),flush=True)
if __name__=='__main__':main()
