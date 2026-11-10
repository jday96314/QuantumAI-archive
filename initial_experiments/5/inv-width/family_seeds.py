"""Independent algorithm seeds: no incumbent arithmetic, model, or tables imported.

Each table is fitted against its own row-wise integer rounding pipeline.
Resources here are explicit analytic circuit schedules (labelled estimates).
"""
import argparse,json,math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from invsqrt.verify_submission_invsqrt import rotation_t_count
MAX_X=3133452

def coordinates(xs,frac):
 k=21-(np.frexp(xs.astype(float))[1]-1);u=xs<<k
 addr=((u>>(21-frac))&((1<<frac)-1))+((k&1)<<frac)
 z=(u&((1<<(21-frac))-1))-(1<<(20-frac))
 return k,u,addr,z

def row_product(a,z,exponent,mode):
 """Signed coordinate, row-rounded products: exact specified integer map."""
 L=21-int(mode.split(':')[1]);kind=mode.split(':')[0]
 out=np.zeros(np.broadcast_shapes(np.shape(a),np.shape(z)),dtype=np.int64)
 if kind=='serial':
  unsigned=z+(1<<(L-1))
  for j in range(L):
   e=exponent+j
   value=a<<e if e>=0 else (a>>-e)+((a>>(-e-1))&1)
   out+=(((unsigned>>j)&1)-(1 if j==L-1 else 0))*value
 elif kind=='booth':
  unsigned=z&((1<<L)-1)
  for j in range((L+1)//2):
   bit=lambda i:0 if i<0 else (unsigned>>min(i,L-1))&1
   digit=bit(2*j-1)+bit(2*j)-2*bit(2*j+1)
   for mag in [1,2]:
    e=exponent+2*j+(mag==2);value=a<<e if e>=0 else (a>>-e)+((a>>(-e-1))&1)
    out+=np.where(abs(digit)==mag,np.sign(digit)*value,0)
 else:raise ValueError(kind)
 return out

def real_coefficients(frac,widths):
 h=2.**-frac;rows=[]
 for parity in range(2):
  scale=1/math.sqrt(2) if parity==0 else 1
  for i in range(1<<frac):
   center=1+(i+.5)*h
   poly=np.polynomial.Chebyshev.interpolate(lambda t:scale/np.sqrt(center+t),len(widths)-1,domain=[-h/2,h/2]).convert(kind=np.polynomial.Polynomial)
   rows.append([poly.coef[j]*(-1)**j*2.**widths[j] for j in range(len(widths))])
 return np.array(rows)

def eval_pipeline(table,z,addr,widths,frac,kind):
 vals=table[addr].T.copy()
 if kind in ('direct','direct_serial'):
  # Direct quadratic: independently rounded slope*z and z^2*curvature.
  square=z*z
  p1=row_product(vals[1],z,widths[0]-widths[1]-21,f'{"serial" if kind=="direct_serial" else "booth"}:{frac}')
  e=widths[0]-widths[2]-42
  p2=vals[2]*square
  p2=(p2+(1<<(-e-1)))>>-e if e<0 else p2<<e
  return vals[0]-p1+p2
 for j in reversed(range(len(widths)-1)):
  vals[j]-=row_product(vals[j+1],z,widths[j]-widths[j+1]-21,f'{"serial" if kind in ("serial","linear_serial") else "booth"}:{frac}')
 return vals[0]

def fit(frac,widths,kind):
 coeff=real_coefficients(frac,widths);table=np.rint(coeff).astype(np.int64)
 # Coordinate descent on quantized high-order coefficients; refit the
 # constant to the midpoint of the *integer-pipeline* error extrema.
 L=21-frac;zz=np.linspace(-(1<<(L-1)),(1<<(L-1))-1,1025).astype(np.int64)
 for a in range(len(table)):
  parity=a>>frac;center=1+((a&((1<<frac)-1))+.5)/2**frac
  desired=(1 if parity else 1/math.sqrt(2))/np.sqrt(center+zz/2.**21)*2.**widths[0]
  aa=np.full_like(zz,a)
  def candidate(row):
   t=table[a].copy();table[a]=row
   out=eval_pipeline(table,zz,aa,widths,frac,kind);table[a]=t
   error=out-desired;delta=round(-(error.max()+error.min())/2)
   return float(np.max(np.abs(error+delta))),delta
  for _ in range(2):
   for j in reversed(range(1,len(widths))):
    best=None
    for d in range(-3,4):
     row=table[a].copy();row[j]+=d
     val,dc=candidate(row);row[0]+=dc
     if best is None or val<best[0]:best=(val,row)
    table[a]=best[1]
   _,dc=candidate(table[a]);table[a,0]+=dc
 return table

def model_check(table,frac,widths,kind):
 max_ratio=0.;worst_x=0;min_margin=1.;signed_low=1.;signed_high=-1.;ranges=[1<<60,-(1<<60)]
 for start in range(1,MAX_X+1,131072):
  xs=np.arange(start,min(MAX_X+1,start+131072),dtype=np.int64);k,u,addr,z=coordinates(xs,frac)
  raw=eval_pipeline(table,z,addr,widths,frac,kind);turns=np.ldexp(raw.astype(float),(k//2)-widths[0]-10)
  target=1/np.sqrt(xs);error=turns-target;tol=target*2.**-24;ratio=np.abs(error)/tol
  idx=ratio.argmax()
  if ratio[idx]>max_ratio:max_ratio=float(ratio[idx]);worst_x=int(xs[idx])
  min_margin=min(min_margin,float(np.min(tol-np.abs(error))));signed_low=min(signed_low,float(np.min(error/tol)));signed_high=max(signed_high,float(np.max(error/tol)))
  ranges=[min(ranges[0],int(raw.min())),max(ranges[1],int(raw.max()))]
 return dict(exhaustive_inputs=MAX_X,max_error_over_tolerance=max_ratio,worst_x=worst_x,min_margin_turns=min_margin,signed_ratio_range=[signed_low,signed_high],output_range=ranges,integer_model_pass=max_ratio<1)

def resource_estimate(table,frac,widths,kind,check):
 # Independent schedules use one carry bit for serial controlled Cuccaro
 # adders; lookup uses unary tree, measured replay deletion (no incumbent).
 sizes=[int(v).bit_length() for v in table.max(axis=0)]
 address=frac+1;norm=210;look=(1<<address)-2
 # Erasure via split-address monomials uses fixed bank cost for all words.
 erasure=sum((1<<n)-n-1 for n in [address//2,address-address//2])
 L=21-frac
 if kind=='direct':
  # Store square, slope, curvature in separate phases. A full schoolbook
  # square is rounded only once: exact multiply into 2L bits, then c2.
  linear_rows=L;square_terms=L*(L+1)//2
  linear_cost=2*sum(3*max(1,min(sizes[1],widths[0]-widths[1]-21+j+sizes[1])) for j in range(L))
  square_cost=2*square_terms*max(1,L//2)+2*L*sizes[2]*3
  arithmetic=linear_cost+square_cost
  Q=max(27+sizes[0]+address-1,27+sizes[1]+19+2,27+2*L+sizes[2]+14+2)
  shift_cost=2*(4*(19+15)-15)+2*(4*(27+15)-15)+2*(4*(14+15)-15)
  rotations=34+29+24;loads=3*look;erase=3*erasure
 elif kind=='serial':
  # Bit-serial Horner with full accumulator storage and one clean carry.
  arithmetic=0
  for j in range(len(widths)-1):
   for bit in range(L):
    e=widths[j]-widths[j+1]-21+bit
    m=max(0,min(sizes[j+1]+e,sizes[j])-max(0,e))
    if m:arithmetic+=2*(3*m+2*max(0,sizes[j]-max(0,e)-m))
  Q=27+sum(sizes)+2
  shift_cost=217;rotations=34;loads=look;erase=erasure
 else:
  # Aligned piecewise linear table; adaptive coarse high intervals use
  # identical leaves coalesced when their integer coefficients agree.
  arithmetic=0
  for bit in range(L):
   e=widths[0]-widths[1]-21+bit;m=max(0,min(sizes[1]+e,sizes[0])-max(0,e))
   if m:arithmetic+=2*(3*m+2*max(0,sizes[0]-max(0,e)-m))
  Q=27+sum(sizes)+2
  shift_cost=217;rotations=34;loads=look;erase=erasure
 T=norm+loads+erase+arithmetic+shift_cost
 margin=check['min_margin_turns'];eps=math.pi*margin/rotations if margin>0 else None
 # A failing approximation has no finite accepted objective; quote the
 # resource objective at the incumbent synthesis budget as an obstruction.
 synthesis=(rotations*rotation_t_count(eps)+1)/4 if eps else rotations*7.38+.25
 return dict(status='ESTIMATE',Q=Q,mean_Toffolis_estimate=T,rotations=rotations,rotation_equivalents_estimate=synthesis,toffoli_equivalents_estimate=T+synthesis,objective_estimate=Q*(T+synthesis),accepted_objective=None,stages_estimate=dict(normalization=norm,lookup=loads,erasure=erase,arithmetic=arithmetic,output_shift=shift_cost),coefficient_bits=sizes,rotation_eps_estimate=eps)

def main():
 rows=[]
 for name,frac,widths,kind in [('direct_quadratic',7,(29,21,12),'direct'),('bit_serial_horner',7,(29,21,12),'serial'),('aligned_table_linear',11,(29,19),'linear')]:
  table=fit(frac,widths,kind);np.save('parameters/'+name+'.npy',table)
  check=model_check(table,frac,widths,kind);resource=resource_estimate(table,frac,widths,kind,check)
  row=dict(family=name,frac=frac,widths=widths,kind=kind,model=check,resources=resource);rows.append(row)
  print(json.dumps(row),flush=True);Path('results/family_first_checkpoint.json').write_text(json.dumps(rows,indent=2)+'\n')
if __name__=='__main__':main()
