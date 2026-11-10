"""28-AND GF(4096) cubing by seven GF(16) cubic evaluations.

Write x=a+b*t+c*t^2 in the established cubic tower. A degree-two
polynomial's cube has degree six, so seven evaluations/interpolation recover
it. GF(16) cubing uses the four-AND cover from round-four cubic_trace_field.
All coordinate transforms are binary linear and all identities are exhaustive.
"""
import json
from pathlib import Path
import kickmix as km
from import_published import load
from published_dqi import isomorphism
from cubic_trace_field import data,cube_terms,inverse_columns,linear
m=load();host=km.GF2Field(12);small=km.GF2Field(4);z=data();fw,bk,tr=isomorphism(m)
columns=tuple(tr(v,fw) for v in z['columns']);back=inverse_columns(columns)
# Invert the degree-six Vandermonde matrix in GF(16).
points=list(range(7));mat=[[small.pow(x,j) for j in range(7)]+[int(i==j) for j in range(7)] for i,x in enumerate(points)]
for j in range(7):
 row=next(i for i in range(j,7) if mat[i][j]);mat[j],mat[row]=mat[row],mat[j]
 c=small.pow(mat[j][j],14);mat[j]=[small.mul(c,v) for v in mat[j]]
 for i in range(7):
  if i!=j:
   c=mat[i][j];mat[i]=[v^small.mul(c,w) for v,w in zip(mat[i],mat[j])]
I=[row[7:] for row in mat];ct,lc,_=cube_terms();terms=[];lincols=[0]*12
for at,p in enumerate(points):
 evalcols=[]
 for nativebit in range(12):
  v=back[nativebit];a,b,c=[v>>(4*i)&15 for i in range(3)]
  evalcols.append(a^small.mul(p,b)^small.mul(small.mul(p,p),c))
 # Li(t) is a native-field element; products with embedded GF16 outputs
 # implement interpolation and reduction modulo the tower polynomial.
 value=0;t=columns[4]
 for degree in range(7):
  value^=m.gmul(tr(linear(z['columns'][:4],I[degree][at]),fw),m.gpow(t,degree))
 outcols=[m.gmul(columns[j],value) for j in range(4)]
 for a,b,o in ct:
  l=sum(((a&v).bit_count()%2)<<i for i,v in enumerate(evalcols));r=sum(((b&v).bit_count()%2)<<i for i,v in enumerate(evalcols))
  out=linear(outcols,o);terms.append((l,r,out))
 for i,e in enumerate(evalcols):lincols[i]^=linear(outcols,linear(lc,e))
assert len(terms)==28
for x in range(4096):
 y=linear(lincols,x)
 for l,r,o in terms:
  if (l&x).bit_count()%2 and (r&x).bit_count()%2:y^=o
 assert y==m.gpow(x,3),(x,y,m.gpow(x,3))
res=dict(status='ALGEBRA_PASS',ANDs=len(terms),field_modulus=(1<<12)^m.gmul(2048,2),terms=terms,linear_columns=lincols,inputs_checked=4096,scope='exact field cube cover; quantum emitter not yet tested')
Path('parameters/cube28.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({k:v for k,v in res.items() if k not in ['terms','linear_columns']}),flush=True)
