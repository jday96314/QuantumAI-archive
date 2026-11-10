"""Classical probe of the q-form BM template and Forney Wronskian relation."""
import json,random
from import_published import load
p=load();mul=p.gmul;inv=p.ginv

def evaluate(a,x):
 y=0
 for c in reversed(a):y=mul(y,x)^c
 return y

def model(S):
 sig=[1]+[0]*35;bank=[1]+[0]*36;d=0
 for r in range(70):
  bank=[bank[-1]]+bank[:-1]
  D=0
  for i in range(p.IMAXV[r]+1):D^=mul(sig[i],S[r-i])
  ctl=D!=0 and d>=0 and r!=69;q=mul(D,sig[0]);m=inv(q) if ctl else q
  if ctl:
   for k in range(p.KMSV[r]+1):
    if r==1 and k==1:continue
    sig[k],bank[k]=bank[k],sig[k]
  for k in range(1 if r==69 else 0,p.KMMV[r]+1):sig[k]^=mul(m,bank[k])
  if r<69:d=-d if ctl else d+1
 return sig,bank,d

if __name__=='__main__':
 rng=random.Random(23891);records=[]
 for weight in [1,2,15,30,31,34,35]:
  positions=rng.sample(range(4095),weight);values=[rng.randrange(1,4096) for _ in positions];S=[]
  for k in range(1,71):
   s=0
   for pos,val in zip(positions,values):s^=mul(val,p.EXP[k*pos%4095])
   S.append(s)
  sig,bank,d=model(S);degree=max(i for i,c in enumerate(sig) if c);bdegree=max(i for i,c in enumerate(bank) if c)
  deff=d if d else -1
  rotation=deff%37
  bank=bank[rotation:]+bank[:rotation]
  rows=[]
  for pos,val in zip(positions,values):
   x=p.EXP[-pos%4095];root=evaluate(sig,x);deriv=evaluate([sig[i+1] if i%2==0 else 0 for i in range(35)],x)
   t=evaluate(bank,x);product=mul(mul(t,deriv),val)
   rows.append(dict(x=x,root=root,scaled_product=product))
  powers=[]
  for exponent in range(-140,141):
   scaled=[mul(row['scaled_product'],p.EXP[-exponent*p.LOG[row['x']]%4095]) for row in rows]
   if len(set(scaled))==1:powers.append((exponent,scaled[0]))
  fixed=[mul(row['scaled_product'],p.EXP[-2*weight*p.LOG[row['x']]%4095]) for row in rows]
  records.append(dict(weight=weight,sigma_degree=degree,bank_degree=bdegree,d=d,rotated_degree=max(i for i,c in enumerate(bank) if c),bank_constant=bank[0],powers=powers[:10],root_errors=sum(row['root']!=0 for row in rows),identity_pass=all(c==1 for c in fixed)))
 print(json.dumps(records,indent=2))
