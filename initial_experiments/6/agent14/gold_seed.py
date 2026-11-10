import json
from pathlib import Path
from import_published import load
m=load()
for q in [2,16,256]:
 columns=[m.gpow(1<<i,q) for i in range(12)];terms=[]
 for a,b,o in m.KT:
  l=sum(1<<i for i in a);rm=sum(1<<i for i in b)
  r=sum(((columns[i]&rm).bit_count()%2)<<i for i in range(12));out=sum(1<<i for i in o)
  terms.append((l,r,out))
 Path(f'parameters/gold{q}_seed.txt').write_text(str(len(terms))+'\n'+''.join(f'{l} {r} {o}\n' for l,r,o in terms))
