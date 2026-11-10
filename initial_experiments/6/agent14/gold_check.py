import json
from pathlib import Path
from import_published import load
m=load()
def get(q):
 data=list(map(int,Path(f'parameters/gold{q}_best.txt').read_text().split()));n=data[0];t=[tuple(data[i:i+3]) for i in range(1,len(data),3)];assert len(t)==n
 def nonlinear(x):
  y=0
  for a,b,o in t:
   if (a&x).bit_count()%2 and (b&x).bit_count()%2:y^=o
  return y
 linear=[m.gpow(1<<i,q+1)^nonlinear(1<<i) for i in range(12)]
 for x in range(4096):
  y=nonlinear(x)
  for i,c in enumerate(linear):
   if x>>i&1:y^=c
  assert y==m.gpow(x,q+1),(q,x,y,m.gpow(x,q+1))
 return t,linear
if __name__=='__main__':
 t,l=get(2);r=dict(status='PASS',power=3,ANDs=len(t),linear_columns=l,inputs=4096,scope='exhaustive algebraic quadratic-cover verification; no circuit yet');Path('results/gold2_check.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
