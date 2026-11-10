"""Original dyadic nonuniform linear fit: fine low / coarse high mantissas."""
import json,math
from pathlib import Path
import numpy as np
from family_seeds import row_product,MAX_X
from invsqrt.verify_submission_invsqrt import rotation_t_count

W0,W1=29,19

def fit_leaf(parity,index,step):
 frac=10-(step==2);h=step/1024;center=1+(index+step/2)/1024;scale=1 if parity else 1/math.sqrt(2);L=21-frac
 poly=np.polynomial.Chebyshev.interpolate(lambda z:scale/np.sqrt(center+z),1,domain=[-h/2,h/2]).convert(kind=np.polynomial.Polynomial)
 base=np.rint([poly.coef[0]*2**W0,-poly.coef[1]*2**W1]).astype(np.int64);zz=np.linspace(-(1<<(L-1)),(1<<(L-1))-1,1025).astype(np.int64);desired=scale/np.sqrt(center+zz/2.**21)*2**W0
 best=None
 for d in range(-4,5):
  slope=base[1]+d;pred=base[0]-row_product(slope,zz,W0-W1-21,f'serial:{frac}');err=pred-desired;dc=round(-(err.max()+err.min())/2);val=np.max(np.abs(err+dc))
  if best is None or val<best[0]:best=(val,[base[0]+dc,slope])
 return best[1]

def run(threshold):
 table=np.zeros((2048,2),dtype=np.int64);centers=np.zeros(2048,dtype=np.int64);fracs=np.zeros(2048,dtype=np.int64);leaves=0
 for parity in range(2):
  i=0
  while i<1024:
   step=2 if i>=threshold else 1;coef=fit_leaf(parity,i,step);table[parity*1024+i:parity*1024+i+step]=coef;centers[parity*1024+i:parity*1024+i+step]=(1<<21)+(i*2048)+step*1024;fracs[parity*1024+i:parity*1024+i+step]=10-(step==2);i+=step;leaves+=1
 maxratio=0.;worst=0;minmargin=1.
 for start in range(1,MAX_X+1,131072):
  xs=np.arange(start,min(MAX_X+1,start+131072));k=21-(np.frexp(xs.astype(float))[1]-1);u=xs<<k;addr=((u>>11)&1023)+((k&1)<<10);z=u-centers[addr];co=table[addr];raw=co[:,0].copy()
  for f in [9,10]:
   mask=fracs[addr]==f;raw[mask]-=row_product(co[mask,1],z[mask],W0-W1-21,f'serial:{f}')
  out=np.ldexp(raw.astype(float),k//2-W0-10);target=1/np.sqrt(xs);tol=target*2.**-24;ratio=np.abs(out-target)/tol;ii=ratio.argmax()
  if ratio[ii]>maxratio:maxratio=float(ratio[ii]);worst=int(xs[ii])
  minmargin=min(minmargin,float(np.min(tol-np.abs(out-target))))
 rot=36;synthesis=rot*rotation_t_count(math.pi*minmargin/rot)/4+.25 if minmargin>0 else None
 # Give the nonuniform design FREE interval selection and optimistic 70Q.
 # Unary lookup: leaves-2. Arithmetic estimate from the earlier 10-bit
 # independent linear schedule, unchanged despite the coarse 12-bit rows.
 T=(leaves-2)+210+114+1390+217
 row=dict(threshold_mantissa=1+threshold/1024,leaves=leaves,exhaustive_inputs=MAX_X,max_error_over_tolerance=maxratio,worst_x=worst,min_margin_turns=minmargin,model_pass=maxratio<1,Q_estimate=70,Toffolis_estimate=T,rotations=rot,rotation_equivalents_estimate=synthesis,objective_estimate=70*(T+synthesis) if synthesis else None,scope='optimistic schedule estimate; interval-selection gates omitted; no native candidate')
 np.savez('parameters/adaptive_linear_'+str(threshold)+'.npz',coefficients=table,centers=centers,fracs=fracs);return row
if __name__=='__main__':
 rows=[]
 for threshold in [896,960]:
  row=run(threshold);rows.append(row);print(json.dumps(row),flush=True);Path('results/adaptive_linear.json').write_text(json.dumps(rows,indent=2)+'\n')
