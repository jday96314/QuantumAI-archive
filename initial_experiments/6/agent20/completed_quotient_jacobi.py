"""Quotient Jacobi synthesis with explicit completion metadata.

The regular prefix has nonzero a_i; a_i=b_i=0 after degree d. Record a_i=0,
convert to unnormalized quotient heads, arm only the first zero head as x,
and choose terminal 0 for completion or 1 for full degree. The totalized
Euclid inverse produces (t^(n-d) F, t^(n-d) G), so the extra roots are at zero.
Metadata stays live until literal inversion, which is valid for arbitrary
heads as well. Interior degree gaps are not decoded by this module.
"""
import quotient_jacobi_synthesis as qj


def synthesize(m,E,W,c,flags,inverse=False):
 n=len(W)//2;F,G=qj.layout(W,c)
 def record():
  for a,f in zip(W[::2],flags):
   with m.eq_test(E,a,0) as zero:E.cx(zero,f)
 def terminal():E.x(c[0]);E.cx(flags[-1],c[0])
 def arm():
  for i,a in enumerate(W[::2]):
   if i:
    E.x(flags[i-1]);E.b.ccx(flags[i],flags[i-1],a[0]);E.x(flags[i-1])
   else:E.cx(flags[0],a[0])
 if inverse:
  qj.poly.extract(m,E,F,G);arm();terminal();qj.convert(m,E,W,True);record()
 else:
  record();qj.convert(m,E,W);terminal();arm();qj.poly.extract(m,E,F,G,True)
 return F,G


def model(m,words):
 flags=[not a for a in words[::2]];prev=1;heads=[]
 for i,(a,b) in enumerate(zip(words[::2],words[1::2])):
  k=m.ginv(m.gmul(int(a),prev or 1));q0=m.gmul(int(b),k or 1)
  heads.append((q0,k^(int(flags[i] and (i==0 or not flags[i-1])))));prev=k
 F,G=[int(not flags[-1])],[]
 for q0,q1 in reversed(heads):
  nn=len(F);out=[0]*(nn+1)
  for i in range(nn-1):out[i]=(G[i] if i<len(G) else 0)^m.gmul(q0,F[i])^(m.gmul(q1,F[i-1]) if i else 0)
  out[-2]=m.gmul(q0,F[-1] or 1)^(m.gmul(q1,F[-2]) if nn>1 else 0)
  out[-1]=m.gmul(q1,F[-1] or 1)
  F,G=out,F
 return F,G,[int(x) for x in flags]
