"""A 42-product GF16 cubic interpolation multiplier, distinct from GF4 CRT."""
from pathlib import Path
import json
import kickmix as km
from extension_inverse import setup
from crt_multiplier import small_terms
m=setup()[0];g=km.GF2Field(4)
t=next(t for t in range(2,4096) if m.gpow(t,4)^t^1==0)
def embed(v):
 ans=0
 for j in range(4):
  if v>>j&1:ans^=m.gpow(t,j)
 return ans
def invert_columns(cols):
 basis={}
 for i,v in enumerate(cols):
  e=1<<i
  while v:
   b=v.bit_length()-1
   if b not in basis:basis[b]=(v,e);break
   u,f=basis[b];v^=u;e^=f
 if len(basis)!=len(cols):return None
 rows=[]
 for j in range(len(cols)):
  v=1<<j;e=0
  while v:b=v.bit_length()-1;u,f=basis[b];v^=u;e^=f
  rows.append(e)
 return rows
for r in range(2,4096):
 cols=[m.gmul(m.gpow(r,j),embed(1<<k)) for j in range(3) for k in range(4)];back=invert_columns(cols)
 if back is not None:break
points=list(range(5))
V=[[g.pow(x,j) for j in range(5)] for x in points];mat=[a+[int(i==j) for j in range(5)] for i,a in enumerate(V)]
for i in range(5):
 j=next(j for j in range(i,7) if mat[j][i]);mat[i],mat[j]=mat[j],mat[i];iv=g.pow(mat[i][i],14);mat[i]=[g.mul(v,iv) for v in mat[i]]
 for j in range(5):
  if j!=i:
   v=mat[j][i];mat[j]=[a^g.mul(v,b) for a,b in zip(mat[j],mat[i])]
I=[a[5:] for a in mat];terms=[]
for col,x in enumerate(points):
 evalcols=[]
 for v in back:
  e=0
  for j in range(3):e^=g.mul((v>>(4*j))&15,g.pow(x,j))
  evalcols.append(e)
 for l,h,o in small_terms(g.modulus):
  left=sum(((v&l).bit_count()%2)<<j for j,v in enumerate(evalcols));right=sum(((v&h).bit_count()%2)<<j for j,v in enumerate(evalcols));out=0
  for d in range(5):out^=m.gmul(m.gpow(r,d),embed(g.mul(I[d][col],o)))
  terms.append((left,right,out))
assert len(terms)==45
for a in range(4096):
 for b in [0,1,2,5,137,2047,4095]:
  got=0
  for l,h,o in terms:
   if (a&l).bit_count()%2 and (b&h).bit_count()%2:got^=o
  assert got==m.gmul(a,b),(a,b,got,m.gmul(a,b))
Path('parameters/cubic_mul_seed.txt').write_text(str(len(terms))+'\n'+''.join(f'{l} {h} {o}\n' for l,h,o in terms));Path('parameters/cubic_mul_quadratic_seed.txt').write_text(str(len(terms))+'\n'+''.join(f'{l} {h<<12} {o}\n' for l,h,o in terms));Path('results/cubic_mul_seed.json').write_text(json.dumps(dict(status='ALGEBRA_PASS',r=r,t=t,ANDs=45,cases=4096*7),indent=2)+'\n');print(dict(r=r,t=t,ANDs=45,cases=4096*7))
