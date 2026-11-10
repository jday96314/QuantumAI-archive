"""Restore the selected affine offset modulo2^26, preserving fixed top bit.

The frozen all-input mantissa model is in [2^26,2^27). Both the pre-restored
word and restored word have bit26 equal to one. Restore only the low26 bits.
The constant 3*2^23 then needs one AND, with its carry erased by measurement.
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
  a,c,d=y[23:26];scratch=head[0]
  if not inv:b.x(a);b.cx(a,c)
  b.init_and(a,c,scratch);b.x(scratch);b.cx(scratch,d);b.x(scratch);b.del_and(a,c,scratch)
  if inv:b.cx(a,c);b.x(a)
 def parity(inv):b.cx(flags[0],y[25])
 def coarse(inv):add(b,x[14:21],y[18:26],subtract=not inv,carry_workspace=head,qubit_limit=limit)
 steps=[constant,parity,coarse]
 for step in (reversed(steps) if inverse else steps):step(inverse)


def build():
 previous=affine.restore_affine;affine.restore_affine=restore_affine
 try:return affine.build(limit=90,keep=4,delay=2,product=17,deferred=3)
 finally:affine.restore_affine=previous


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',default='candidates/fixed_top_affine.kmb');a=p.parse_args();path=Path(a.output);c=compact(build(),path);c.to_file(str(path),format='kmb');report=dict(status='BUILT',scope='restoration modulo2^26 with fixed normalized top bit; requires full native verification',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,output=str(path));Path('results/fixed_top_affine_build.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
