"""Build the gap-aware decoder with trimmed, split-optimized products.

Read-only dependencies retain upstream Apache-2.0 attribution in round 3's
import_published.py and round 5's dqi-gap/RESULTS.md. Writes are local.
"""
import argparse,hashlib,json,time
from pathlib import Path
parser=argparse.ArgumentParser();parser.add_argument('--gold',action='store_true');args=parser.parse_args()
from triangular_dp import DPArithmetic
if args.gold:
 from gold_arithmetic import GoldArithmetic
 Arithmetic=GoldArithmetic
else:Arithmetic=DPArithmetic
from mobius_series import BoundedArithmetic
# Affect direct local-head constructors as well as the outer decoder.
for name in ['sparse_middle','offset_middle','triangular','triangular_suffix']+(['product','frobenius','cube','cube_terms'] if args.gold else []):
 setattr(BoundedArithmetic,name,getattr(Arithmetic,name))
import hybrid_decoder as decoder
from blocked_monic_euclid import schedule

counts=json.loads(Path('../../5/dqi-gap/parameters/best_schedule.json').read_text())['counts']
def frames(S,width):
 bp=[];b=[S[-2]];c=S[:-2];p=S[-1];k=1;stages=[]
 for count in counts:
  stages.append((k,count,bp,b,c,p));_,(bp,b)=schedule(bp,b,list(reversed(c))[:2*count]);c=c[:-2*count];k+=count
 assert k==35 and not c
 return stages,(bp,b,c,p)
decoder.ROOT=Path.cwd();decoder.frames=frames;decoder.build.slot_counts=counts
tag='gap_gold' if args.gold else 'gap_dp'
out=Path(f'candidates/{tag}_n70_m4095.kmb')
c=decoder.build(70,3,100,4095,out,linear_padding=True,virtual_zero=True,virtual_lead=True,low_inverse=True,lazy_flags=True)
h=hashlib.sha256(out.read_bytes()).hexdigest()
res=dict(status='BUILT',scope='requires full default native verification',Q=c.num_qubits,T=c.operation_statistics().num_toffoli_like,sha256=h,counts=counts)
Path(f'results/{tag}_manifest.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res),flush=True)
