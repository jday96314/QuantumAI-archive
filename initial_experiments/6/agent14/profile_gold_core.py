import json
from pathlib import Path
from gold_arithmetic import GoldArithmetic
from mobius_series import BoundedArithmetic
for name in ['sparse_middle','offset_middle','triangular','triangular_suffix','product','frobenius','cube','cube_terms']:
 setattr(BoundedArithmetic,name,getattr(GoldArithmetic,name))
import hybrid_decoder as decoder
from profile_hybrid import module
from blocked_monic_euclid import schedule
counts=json.loads(Path('../../5/dqi-gap/parameters/best_schedule.json').read_text())['counts']
def frames(S,width):
 bp=[];b=[S[-2]];c=S[:-2];p=S[-1];k=1;stages=[]
 for count in counts:
  stages.append((k,count,bp,b,c,p));_,(bp,b)=schedule(bp,b,list(reversed(c))[:2*count]);c=c[:-2*count];k+=count
 return stages,(bp,b,c,p)
decoder.ROOT=Path.cwd();decoder.frames=frames;decoder.load=module;decoder.build.profile_only=True
r=decoder.build(70,3,100,4095,Path('candidates/gap_gold_profile.kmb'),True,True,True,True,True)
Path('results/gold_core_profile.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='rows'}),flush=True)
