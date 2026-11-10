"""Generate an exact14-product GF4 tensor for GF4096 via CRT.

Uses the existing independently checked projective-line polynomial helpers;
keeps GF4 operations whole for a search with nonbinary projective scalings.
"""
from pathlib import Path
import kickmix as km
import gf4_multiplier as g
from tower_inverse import map_ops


def seed():
 f=km.GF2Field(12);root=next(x for x in range(2,4096) if f.mul(x,x)^x^1==0);embed=[0,1,root,root^1];theta=f.primitive_element
 ops=map_ops([f.mul(embed[1<<j],f.pow(theta,i)) for i in range(6) for j in range(2)])
 def evaluate(poly):
  out=0
  for i in reversed(range(g.degree(poly)+1)):out=f.mul(out,theta)^embed[(poly>>(2*i))&3]
  return g.map_value(out,ops,True)
 quads=[16|(p<<2)|q for p in range(4) for q in range(1,4) if all(g.remainder(16|(p<<2)|q,4|a) for a in range(4))]
 mods=[4|a for a in range(4)]+quads[:3];P=1
 for m in mods:P=g.multiply(P,m)
 terms=[]
 for m in mods:
  q=g.quotient(P,m);E=g.multiply(q,g.invert(g.remainder(q,m),m));values=[g.remainder(1<<(2*i),m) for i in range(6)];a=sum((v&3)<<(2*i) for i,v in enumerate(values));b=sum(((v>>2)&3)<<(2*i) for i,v in enumerate(values))
  def image(poly):return evaluate(g.remainder(poly,P))
  if g.degree(m)==1:terms.append((a,a,image(E)))
  else:
   p=(m>>2)&3;q=m&3
   terms.extend([(a,a,image(g.multiply(E,5))),(b,b,image(g.multiply(E,q|((1^p)<<2)))),(a^b,a^b,image(g.multiply(E,4)))])
 terms.append((1<<10,1<<10,evaluate(P)));assert len(terms)==14
 for i in range(6):
  for j in range(6):
   out=0
   for a,b,c in terms:out^=g.scalar(c,g.F4.mul((a>>(2*i))&3,(b>>(2*j))&3))
   assert out==g.map_value(f.mul(f.pow(theta,i),f.pow(theta,j)),ops,True)
 return terms


if __name__=='__main__':
 terms=seed();Path('parameters/gf4_rank14_initial.txt').write_text(str(len(terms))+'\n'+''.join(f'{a} {b} {c}\n' for a,b,c in terms));print('PASS:14GF4 products, all36 field-basis pairs')
