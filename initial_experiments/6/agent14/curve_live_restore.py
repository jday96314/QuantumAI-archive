"""Exact curvature replay from the restored affine C0 and coarse address."""
import json
from pathlib import Path
import kickmix as km


def apply(b,c0,address,slope,high,keep=0):
 assert keep==0
 z=json.loads(Path('parameters/curve_live_cut13.json').read_text());words=list(c0[:21])+list(address[1:])+list(slope[13:]);assert len(words)==z['width'];outs=list(high)
 def linear(mask,target):
  for i,q in enumerate(words):
   if mask>>i&1:b.cx(q,target)
  if mask>>len(words)&1:b.x(target)
 for q,mask in zip(outs,z['linear']):linear(mask,q)
 def product(l,r,targets):
  constl=l>>len(words)&1;constr=r>>len(words)&1;left={i for i in range(len(words)) if l>>i&1};right={i for i in range(len(words)) if r>>i&1}
  if constl:
   for t in targets:linear(r,t)
  if constr:
   for t in targets:
    for i in left:b.cx(words[i],t)
  if not left or not right:return
  if left==right:
   for t in targets:
    for i in left:b.cx(words[i],t)
   return
  p=min(left);first=sorted(left-{p})
  for i in first:b.cx(words[i],words[p])
  if p in right:right.symmetric_difference_update(left-{p})
  q=min(right-{p});second=sorted(right-{q})
  for i in second:b.cx(words[i],words[q])
  for t in targets[1:]:b.cx(targets[0],t)
  b.ccx(words[p],words[q],targets[0])
  for t in targets[1:]:b.cx(targets[0],t)
  for i in reversed(second):b.cx(words[i],words[q])
  for i in reversed(first):b.cx(words[i],words[p])
 for l,r,o in zip(z['left'],z['right'],z['output']):product(l,r,[q for j,q in enumerate(outs) if o>>j&1])
 # Fine leaf69 has the sole C0/coarse collision; the fitted base uses leaf68.
 correction=(37^38)>>keep
 targets=[q for j,q in enumerate(outs) if correction>>j&1]
 for q in targets[1:]:b.cx(targets[0],q)
 if targets:b.flip_if_equal(address,69,target=targets[0])
 for q in targets[1:]:b.cx(targets[0],q)
