"""Six-moment coordinates using two nonsingular adjacent Padé minors.

The map allows a zero first Euclidean pivot. Its inverse is explicit; it does
not invoke that singular pivot. This is a classical research component only,
without a reversible quantum schedule or a complete decoder resource claim.
"""
import argparse,json,random,time
from pathlib import Path
from compact_euclid import mul,inv,divide,xor,product


def add(*terms):
 value=0
 for v in terms:value^=v
 return value


def square(a):return mul(a,a)


def forward(F):
 e,d,a,b,c,f=F
 D2=square(d)^mul(e,a)
 if not D2:raise ValueError('singular second Pade minor')
 D2inv=inv(D2);u0=mul(mul(d,b)^square(a),D2inv);u1=mul(mul(e,b)^mul(d,a),D2inv)
 p=add(c,mul(a,u0),mul(b,u1))
 if not p:raise ValueError('singular third Pade minor')
 alpha=mul(add(f,mul(d,square(u0)),mul(b,square(u1))),inv(p));z=mul(p,D2inv)
 v0=mul(alpha,u0)^mul(z,d);v1=add(u0,mul(alpha,u1),mul(z,e));v2=u1^alpha
 return [u0,u1,v0,v1,v2,inv(p)]


def backward(C):
 u0,u1,v0,v1,v2,pinv=C
 if not pinv:raise ValueError('zero reciprocal normalizer')
 alpha=v2^u1;r0=v0^mul(alpha,u0);r1=add(v1,u0,mul(alpha,u1))
 z=mul(add(square(r0),mul(square(r1),u0),mul(mul(r0,r1),u1)),pinv)
 if not z:raise ValueError('singular inverse Pade minor')
 zi=inv(z);e=mul(r1,zi);d=mul(r0,zi);a=mul(e,u0)^mul(d,u1);b=mul(d,u0)^mul(a,u1);p=inv(pinv)
 c=add(p,mul(a,u0),mul(b,u1));f=add(mul(alpha,p),mul(d,square(u0)),mul(b,square(u1)))
 return [e,d,a,b,c,f]


def transfer(F,C):
 e,d,a,b,_,_=F;u0,u1,v0,v1,v2,pinv=C
 A2=[d^mul(e,u1),e];B2=[u0,u1,1]
 A3=[add(a,mul(d,v2),mul(e,v1)),d^mul(e,v2),e];B3=[v0,v1,v2,1]
 det=xor(product(A2,B3),product(A3,B2))
 assert det[0]==inv(pinv) and not any(det[1:]),(F,C,det)
 return A2,B2,A3,B3


def main():
 ap=argparse.ArgumentParser();ap.add_argument('--shots',type=int,default=12000);a=ap.parse_args();rng=random.Random(152731);begin=time.time();counts={};rejected={};steps={};checks=0;inverse_checks=0
 for shot in range(a.shots):
  F=[rng.randrange(4096) for _ in range(6)];kind=['random','zero_first','zero_second','zero_first_two'][shot%4]
  if kind=='zero_first':F[0]=0
  if kind=='zero_second':F[1]=0
  if kind=='zero_first_two':F[0]=F[1]=0
  try:C=forward(F)
  except ValueError as ex:rejected[kind]=rejected.get(kind,0)+1;continue
  assert backward(C)==F,(kind,F,C,backward(C));inverse_checks+=1;A2,B2,A3,B3=transfer(F,C)
  P=[0]*6+[1];R=list(reversed(F));bp=[];bb=[1];degree_steps=[];found=[]
  while R and not R[-1]:R.pop()
  while len(bb)<=4 and R:
   q,nxt=divide(P,R);newb=xor(bp,product(q,bb));P,R=R,nxt;bp,bb=bb,newb
   while bb and not bb[-1]:bb.pop()
   lead=inv(bb[-1]);monic=[mul(lead,v) for v in bb];degree_steps.append(len(q)-1)
   if len(monic)==3:found.append(monic)
   if len(monic)==4:
    assert found and found[-1]==B2,(kind,F,C,found,B2)
    assert monic==B3,(kind,F,C,monic,B3)
    break
  else:raise AssertionError((kind,F,C,'missing degree3 cofactor'))
  key=','.join(map(str,degree_steps));steps[key]=steps.get(key,0)+1;counts[kind]=counts.get(kind,0)+1;checks+=1
 result=dict(status='PASS',scope='classical invertible six-moment/Pade-coordinate map, full Euclid comparison, transfer-matrix determinant, including first-pivot zero; no quantum resources or full decoding acceptance',requested_shots=a.shots,checks=checks,inverse_checks=inverse_checks,accepted_by_kind=counts,rejected_by_kind=rejected,euclid_quotient_degree_patterns=steps,seconds=time.time()-begin)
 Path('results/lookahead_six_moment_map.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
