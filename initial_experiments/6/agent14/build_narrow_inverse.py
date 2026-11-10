"""Compile the new reduced-curvature table into the proven inverse schedule."""
import argparse,inspect,json,math
from pathlib import Path
import numpy as np
import kickmix as km
import affine_invsqrt as affine
import binary_centered_invsqrt as polynomial
import fixed_top_affine_invsqrt as fixed
import zero_head_affine_invsqrt as zero
from optimize_workspace import compact
from invsqrt.verify_submission_invsqrt import rotation_t_count

TABLE=np.loadtxt('parameters/narrow_fit_27_18_7.txt',dtype=np.int64,skiprows=1)
def affine_table():
 t=TABLE.copy();idx=np.arange(len(t));parity=idx>>8;coarse=(idx&255)>>1
 t[:,0]+=(1<<25)+(coarse<<18)+(1<<17)-np.array(affine.A)[parity]
 return t

def build(limit,keep,delay):
 src=inspect.getsource(affine.build).replace('(27,18,8)','(27,18,7)')
 ns=dict(affine.__dict__);ns.update(ROOT=Path.cwd(),table=affine_table,restore_affine=fixed.restore_affine)
 exec(src,ns);previous=polynomial.cshift;polynomial.cshift=zero.cshift
 try:return ns['build'](limit=limit,keep=keep,delay=delay,product=17,deferred=3)
 finally:polynomial.cshift=previous

def numerical():
 from invsqrt.verify_submission_invsqrt import ProblemSpecInvSqrt,reference_values
 spec=ProblemSpecInvSqrt.for_n_p(10);previous=polynomial.table_for;polynomial.table_for=lambda *args:TABLE
 worst=0;minmargin=1<<64;worstx=0;where=0
 try:
  for first in range(1,3133453,131072):
   xs=np.arange(first,min(first+131072,3133453),dtype=np.int64)
   p=polynomial.model(xs,7,(27,18,7),rounded=True,lookup_frac=8,booth=True).astype(np.uint64)<<np.uint64(27)
   target,tol=reference_values(spec,xs.astype(np.uint64));used=np.abs((p-target).view(np.int64)).astype(np.uint64)
   ratio=used/tol;j=int(ratio.argmax())
   if float(ratio[j])>worst:worst=float(ratio[j]);worstx=int(xs[j])
   margin=tol.astype(np.int64)-used.astype(np.int64);j=int(margin.argmin())
   if int(margin[j])<minmargin:minmargin=int(margin[j]);where=int(xs[j])
 finally:polynomial.table_for=previous
 res=dict(status='MODEL_PASS' if minmargin>0 else 'MODEL_FAIL',inputs=3133452,worst_ratio=worst,worst_x=worstx,min_margin_units=minmargin,margin_x=where)
 Path('results/narrow_inverse_model.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res),flush=True)
 return res

if __name__=='__main__':
 model=numerical();assert model['status']=='MODEL_PASS';eps=math.pi*model['min_margin_units']/(1<<64)/34;rot=(1+34*rotation_t_count(eps))/4
 rows=[]
 for limit in [86,88,90]:
  for keep in [2,3,4]:
   try:
    c=compact(build(limit,keep,2),Path('candidates/narrow_inverse_temp.kmb'));T=c.operation_statistics().num_toffoli_like
    row=dict(limit=limit,keep=keep,Q=c.num_qubits,T=T,rotation_equivalents=rot,estimated_objective=c.num_qubits*(T-11+rot));rows.append(row)
    if row['estimated_objective']==min(r['estimated_objective'] for r in rows):c.to_file('candidates/narrow_inverse_best.kmb',format='kmb')
    print(json.dumps(row),flush=True)
   except (ValueError,AssertionError) as e:print(json.dumps(dict(limit=limit,keep=keep,error=str(e))),flush=True)
 Path('results/narrow_inverse_builds.json').write_text(json.dumps(rows,indent=2)+'\n')
