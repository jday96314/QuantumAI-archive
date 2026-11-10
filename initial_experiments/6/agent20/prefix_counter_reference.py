"""Independent scalar-field reference for the complete research oracle."""
import check_jacobi
import check_prefix_high_copy as early
import check_high_copy_jacobi as late
import prefix_high_copy_jacobi as schedule

def extract(F,syndrome,pairs=3,through=3,suffix=2):
 class Adapter:
  gmul=staticmethod(F.mul);gpow=staticmethod(F.pow);ginv=staticmethod(lambda a:F.inv(a) if a else 0)
 W=list(syndrome);k=0;early_regular=True
 for b in schedule.schedule(len(W)//2,pairs,through,suffix):
  if b==1:W[2*k:]=check_jacobi.pair(Adapter,W[2*k:])
  else:
   module=early if k<through else late;module.PAIRS=b
   if k<through and len(W[2*k:])>2*b:
    H=late.head(Adapter,W[2*k:2*(k+b)])
    A,_,_=early.frame_model(Adapter,H)
    early_regular &= bool(A[-1]) and all(H[::2])
   W[2*k:]=module.model(Adapter,W[2*k:])
  k+=b
 return W,early_regular

def synthesize(F,W):
 W=list(W);n=len(W)//2;zeros=sum(a==0 for a in W[::2]);prev=1
 for i in range(n):
  a,b=W[2*i:2*i+2];p=F.mul(a,prev or 1);k=F.inv(p) if p else 0
  W[2*i:2*i+2]=[k,F.mul(b,k or 1)];prev=k
 d=n-zeros
 if zeros and W[2*d]==0:W[2*d]=1
 f,g=[int(zeros==0)],[]
 for q1,q0 in reversed(list(zip(W[::2],W[1::2]))):
  out=[0]*(len(f)+1)
  for i,v in enumerate(g):out[i]=v^F.mul(q0,f[i])^(F.mul(q1,f[i-1]) if i else 0)
  out[-2]=F.mul(q0,f[-1] or 1)^(F.mul(q1,f[-2]) if len(f)>1 else 0)
  out[-1]=F.mul(q1,f[-1] or 1);f,g=out,f
 return f,g

def model(F,syndrome,pairs=3,through=3,suffix=2):
 W,regular=extract(F,syndrome,pairs,through,suffix);f,g=synthesize(F,W)
 return f,g,regular

def phase(F,f,g,masks):
 value=0;derivative=F.derivative(f)
 for j,mask in enumerate(masks):
  a=F.pow(2,j%4095)
  if F.evaluate(f,a)==0:
   den=F.mul(a,F.evaluate(derivative,a));e=F.mul(F.evaluate(g,a),F.inv(den)) if den else 0
   value^=(e&int(mask)).bit_count()%2
 return value
