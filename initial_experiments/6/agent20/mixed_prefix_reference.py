"""Independent totalized scalar reference for mixed prefix/suffix oracles."""
import prefix_counter_reference as prior
import check_high_copy_jacobi as headref
import check_prefix_high_copy as prefixref
import check_jacobi
import check_blocked_cfraction as ops
import prefix_high_copy_jacobi as schedule

def look_block(m,words):
 H,R=list(words[:4]),list(words[4:]);flag=not H[0];n=len(R)
 if flag:
  a,u,v,w=H;e=m.gmul(v,m.ginv(u or 1));h=m.gmul(w,m.ginv(u or 1))^m.gmul(e,e);H=[a,u,h,e]
 else:
  H=headref.head(m,H);a,b,c,d=H;H=[a,b,c^m.gmul(b,d),b^d]
 a,b,h,e=H;c=h^m.gmul(b,e)^m.gmul(b,b);d=e^b;bg=b if flag else 0;A=[1,e,h];D=[1,b^bg];p=(b if flag else m.gmul(a,c or 1)) or 1
 P=[m.gmul(a,m.gpow(b,4)^m.gmul(c,c)^m.gmul(c,m.gmul(e,e)))^m.gmul(bg,m.gpow(e,3)),m.gmul(h,m.gmul(a,m.gpow(b,3)^m.gmul(c,d))^m.gmul(bg,h^m.gmul(e,e)))];shift=2+flag
 if n:
  R=ops.product(m,A,R,n)
  for i,v in enumerate(P[:n]):R[i]^=v
  R=ops.divide(m,ops.product(m,A,R,n),[p,*([0]*(shift-1)),*ops.product(m,D,R,max(n-shift,0))],n)
 return H+R

def synthesize(F,W,through):
 W=list(W);n=len(W)//2;flags=[];zeros=0;prev=1
 for k in range(0,n-1,2):
  a,b,h,e=W[2*k:2*k+4];flag=k>=through and a==0;flags.append(flag)
  if k>=through:
   if flag:c,d=e,h
   else:c,d=h^F.mul(b,e^b),e^b
  else:c,d=h,e
  if flag:
   zeros+=2*(b==0);v=F.mul(b,prev or 1);q=F.inv(v) if v else 0;W[2*k:2*k+4]=[q,F.mul(c,q or 1),0,F.mul(d,q or 1)];prev=q
  else:
   zeros+=(a==0)+(c==0);v=F.mul(a,prev or 1);q=F.inv(v) if v else 0;v=F.mul(c,q or 1);r=F.inv(v) if v else 0;W[2*k:2*k+4]=[q,F.mul(b,q or 1),r,F.mul(d,r or 1)];prev=r
 if n%2:
  a,b=W[-2:];zeros+=a==0;v=F.mul(a,prev or 1);q=F.inv(v) if v else 0;W[-2:]=[q,F.mul(b,q or 1)]
 d=n-zeros
 if zeros and W[2*d]==0:W[2*d]=1
 def one(q0,q1,f,g):
  out=[0]*(len(f)+1)
  for i,v in enumerate(g):out[i]=v^F.mul(q0,f[i])^(F.mul(q1,f[i-1]) if i else 0)
  out[-2]=F.mul(q0,f[-1] or 1)^(F.mul(q1,f[-2]) if len(f)>1 else 0);out[-1]=F.mul(q1,f[-1] or 1)
  return out,f
 f,g=[int(zeros==0)],[]
 if n%2:f,g=one(W[-1],W[-2],f,g)
 for packet in reversed(range(n//2)):
  q2,q1,r1,q0=W[4*packet:4*packet+4]
  if not flags[packet]:
   f,g=one(q0,r1,f,g);f,g=one(q1,q2,f,g)
  else:
   gg=[*f,0];out=[*g,q0,q1,q2];nn=len(gg)
   for i in range(nn-2):
    for j,q in enumerate((q0,q1,q2)):
     if i>=j:out[i]^=F.mul(q,gg[i-j])
   lead=gg[-2] or 1;out[-3]=F.mul(q0,lead)^(F.mul(q1,gg[-3]) if nn>=3 else 0)^(F.mul(q2,gg[-4]) if nn>=4 else 0);out[-2]=F.mul(q1,lead)^(F.mul(q2,gg[-3]) if nn>=3 else 0);out[-1]=F.mul(q2,lead);f,g=out,gg
 return f,g

def model(F,s,pairs=2,through=2,suffix=2):
 class Adapter:
  gmul=staticmethod(F.mul);gpow=staticmethod(F.pow);ginv=staticmethod(lambda a:F.inv(a) if a else 0)
 W=list(s);k=0;regular=True
 for b in schedule.schedule(len(W)//2,pairs,through,suffix):
  if b==1:W[2*k:]=check_jacobi.pair(Adapter,W[2*k:])
  elif k<through:
   prefixref.PAIRS=b;H=headref.head(Adapter,W[2*k:2*(k+b)]);A,_,_=prefixref.frame_model(Adapter,H);regular &= bool(A[-1]) and all(H[::2]);W[2*k:]=prefixref.model(Adapter,W[2*k:])
  else:W[2*k:]=look_block(Adapter,W[2*k:])
  k+=b
 f,g=synthesize(F,W,through);return f,g,regular

phase=prior.phase
