"""Full-EEA diagnosis of the retained law-sample failure; no decoder changes."""
import json
from pathlib import Path
from matrix_transfer import native_reference

mul,inv,_,_=native_reference()
row=json.loads(Path('results/law_1000000_2026102405_g3_failures.jsonl').read_text())
p=[0]*70+[1];r=row['syndromes'].copy();stage=0;gaps=[]
while len(r)>35:
    g=len(p)-len(r)
    if g>1:gaps.append(dict(degree_stage=stage,quotient_degree=g))
    a=p.copy();ri=inv(r[-1])
    while len(a)>=len(r):
        d=len(a)-len(r);v=mul(a[-1],ri)
        for j,x in enumerate(r):a[d+j]^=mul(x,v)
        while a and not a[-1]:a.pop()
    p,r=r,a;stage+=g
assert gaps==[dict(degree_stage=24,quotient_degree=3)]
report=dict(shot=row['shot'],seed=2026102405,gaps=gaps,completed_locator_degree=stage,
            evaluator_degree=len(r)-1,
            diagnosis='degree-three quotient at stage 24 exceeds two-slot serial budget',
            scope='independent full polynomial EEA stopped at the decoding bound')
Path('results/serial_failure_diagnosis.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
