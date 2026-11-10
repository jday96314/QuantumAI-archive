"""Use measured cleanup for barrel-shift moves into known-zero headroom.

A move (s,0)->((1-c)s,cs) uses AND then CNOT. Its inverse restores s by
CNOT and erases cs by X-basis measurement with the AND phase correction.
Only private phase headroom is measured; normalization retains ordinary swaps.
"""
import argparse,json
from pathlib import Path
import kickmix as km
import binary_centered_invsqrt as polynomial
import fixed_top_affine_invsqrt as fixed_top
from optimize_workspace import compact

original_shift=polynomial.cshift


def cshift(b,reg,k,ctrl,occ,inv=False):
 if len(reg)==22:return original_shift(b,reg,k,ctrl,occ,inv)
 assert len(reg)==37
 indices=[i for i in range(min(occ,len(reg))) if i+k<len(reg)]
 for i in (indices if inv else reversed(indices)):
  source,target=reg[i],reg[i+k]
  if i==26 and k==1 and occ==27:
   # The selected all-x model proves this source is the constant one.
   b.cx(ctrl,source);b.cx(ctrl,target)
  elif i+k>=occ:
   if inv:b.cx(target,source);b.del_and(ctrl,source,target)
   else:b.init_and(ctrl,source,target);b.cx(target,source)
  else:b.cswap(ctrl,source,target)


def build():
 previous=polynomial.cshift;polynomial.cshift=cshift
 try:return fixed_top.build()
 finally:polynomial.cshift=previous


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',default='candidates/zero_head_affine.kmb');a=p.parse_args();path=Path(a.output);c=compact(build(),path);c.to_file(str(path),format='kmb');report=dict(status='BUILT',scope='fixed-top affine restoration and measured clean-headroom shifts; requires full native verification',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,output=str(path));Path('results/zero_head_affine_build.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
