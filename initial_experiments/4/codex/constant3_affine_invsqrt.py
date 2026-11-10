"""Combine the selected affine restoration's adjacent constant increments.

Adding 3 at bit23 uses one two-bit addition with a repeated constant-one
source. The adder has sufficient borrowed carry space to use its plain path,
which restores every source before another aliased source is read.
"""
import argparse,json
from pathlib import Path
import kickmix as km
import affine_invsqrt as affine
from short_adder import add
from optimize_workspace import compact


def restore_affine(b,x,y,regs,flags,limit,inverse=False):
 head=y[27:]
 assert len(head)==10 and x[21]==y[26]
 def constant(inv):
  one=head[-1];b.x(one)
  assert len(head[:-1])>=len(y[23:27])-2
  add(b,km.array([one,one]),y[23:27],subtract=inv,carry_workspace=head[:-1],qubit_limit=limit)
  b.x(one)
 def parity(inv):add(b,km.array([flags[0]]),y[25:27],subtract=inv,carry_workspace=head,qubit_limit=limit)
 def coarse(inv):add(b,x[14:21],y[18:27],subtract=not inv,carry_workspace=head,qubit_limit=limit)
 steps=[constant,parity,coarse]
 for step in (reversed(steps) if inverse else steps):step(inverse)


def build():
 previous=affine.restore_affine;affine.restore_affine=restore_affine
 try:return affine.build(limit=90,keep=4,delay=2,product=17,deferred=3)
 finally:affine.restore_affine=previous


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',default='candidates/constant3_affine.kmb');a=p.parse_args();path=Path(a.output);c=compact(build(),path);c.to_file(str(path),format='kmb');report=dict(status='BUILT',scope='combined constant restoration; requires full native verification',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,output=str(path));Path('results/constant3_affine_build.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
