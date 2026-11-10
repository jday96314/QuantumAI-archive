"""All coefficients of independent C++ and native-validated scalar models."""
from pathlib import Path
import hashlib,json,subprocess
import numpy as np
import mixed_prefix_reference as reference
from algebra import Field
HERE=Path(__file__).resolve().parent
F=Field();assert F.mul(2048,2)==0x53
rng=np.random.default_rng(2026102622);fixtures=[];expected=[];counts={}
for n,p,t,trials in ((12,2,2,4),(14,2,2,4),(70,6,18,2)):
 count=0
 for d in range(n//2+1):
  for trial in range(trials):
   pos=rng.choice(4095,d,replace=False);roots=[F.pow(2,int(j)) for j in pos];values=rng.integers(1,4096,d).tolist();s=F.syndrome(roots,values,n)
   f,g,_=reference.model(F,s,p,t,2);fixtures.append(' '.join(map(str,[n,p,t,*s])));expected.append(f+g);count+=1
 for i in range(16):
  s=rng.integers(0,4096,n).tolist()
  if i%2==0:s[::2]=[0]*(n//2)
  f,g,_=reference.model(F,s,p,t,2);fixtures.append(' '.join(map(str,[n,p,t,*s])));expected.append(f+g);count+=1
 counts[n]=count
r=subprocess.run([str(HERE/'mixed_prefix_model'),'--fixtures'],input='\n'.join(fixtures)+'\n',text=True,check=True,capture_output=True)
actual=[list(map(int,line.split())) for line in r.stdout.splitlines()];assert actual==expected,next(((i,a,b) for i,(a,b) in enumerate(zip(actual,expected)) if a!=b),'length mismatch')
sha=lambda p:hashlib.sha256((HERE/p).read_bytes()).hexdigest()
report=dict(status='PASS_INDEPENDENT_MODEL_AGREEMENT',cases=len(fixtures),by_N=counts,seed=2026102622,source_sha256={p:sha(p) for p in ['mixed_prefix_model.cpp','prefix_counter_model.cpp','generalized_quotient_jacobi_model.cpp','mixed_prefix_reference.py','check_mixed_prefix_model.py']},scope='all polynomial coefficients, every completion degree at N12/N14/N70, arbitrary and forced nonunary zero inputs; Python model separately matched native small oracle; no native full-task acceptance')
(HERE/'results/mixed_prefix_model_agreement.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
