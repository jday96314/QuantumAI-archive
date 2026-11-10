"""Full selected decoder with a forced stage-10 gap branch.

Diagnostic only: requires a single degree gap exactly at stage 10 and a
nonzero following normalizer. This is NOT a submission candidate or a global
repair. It tests the new branch through the selected prefix, phase scan and
inverse on the original syndrome contract. A global quantum dispatcher is absent.
"""
import json,hashlib
from pathlib import Path
import compact_product_decoder as decoder
from cubic_inverse_product_kernel import apply as phase_kernel
from gap_block import apply as gap
from blocked_monic_euclid import schedule
from quantum_coprime_padding import permute
regular=decoder.block

def branch(m,E,math,Bp,B,C,pinv,count,inverse=False,**kw):
 if len(B)!=10 or count!=3:return regular(m,E,math,Bp,B,C,pinv,count,inverse,**kw)
 raw=list(reversed(C))[:6];Q=[raw[i] for i in [1,3,2,4,5,0]]
 gp=Bp+Q[:3];gb=B+[Q[3],Q[5],Q[4]]
 _,(np,nb)=schedule(Bp,B,raw);ops=[];permute(ops,gp+gb,np+nb)
 def align(back=False):
  for op in (reversed(ops) if back else ops):getattr(E.b,op[0])(*op[1:])
 if inverse:align(True);gap(m,E,math,Bp,B,C,pinv,True)
 else:gap(m,E,math,Bp,B,C,pinv);align()
 return np,nb,C[:-6],pinv

decoder.block=branch;decoder.phase_kernel=phase_kernel;decoder.ROOT=Path.cwd();out=Path('candidates/gap10_diagnostic.kmb')
c=decoder.build(70,3,22,4095,out,linear_padding=True,virtual_zero=True,virtual_lead=True,low_inverse=True,lazy_flags=True)
manifest=dict(scope=__doc__,sha256=hashlib.sha256(out.read_bytes()).hexdigest(),files={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in ['build_gap10_diagnostic.py','gap_block.py','gap_cofactor.py']})
Path('parameters/gap10_diagnostic_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
