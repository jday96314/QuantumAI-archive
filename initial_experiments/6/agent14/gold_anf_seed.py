from pathlib import Path
from import_published import load
m=load()
for q in [2,16]:
 terms=[]
 for i in range(12):
  for j in range(i):
   o=m.gmul(1<<i,m.gpow(1<<j,q))^m.gmul(1<<j,m.gpow(1<<i,q))
   if o:terms.append((1<<i,1<<j,o))
 Path(f'parameters/gold{q}_anf.txt').write_text(str(len(terms))+'\n'+''.join(f'{a} {b} {o}\n' for a,b,o in terms))
