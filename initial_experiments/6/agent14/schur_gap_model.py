"""History-free scalar-pivot Schur chart: independent classical prototype.

An isolated b=0 uses two reciprocal steps with fixed nodes. Unlike swapping,
the first zero b remains a chart tag in the stored field words. This model
stops completed zero tails explicitly; a reversible full-width circuit and
in-place mixed-degree synthesis are still required.
"""
import json,random,time
from pathlib import Path
from quadratic_schur_model import F,K,node_maps,transform,bmul,einv,norm


def extract(N,syndrome):
 n=N//2;nodes,kap,P,maps=node_maps(n);v=transform(syndrome,maps);g=1;stages=[];exc=None
 for i in range(n):
  if not i:A=B=0
  else:
   pa,pb=stages[-1]['a'],stages[-1]['b']
   if not pb:A=B=0
   elif i>1 and not stages[-2]['b']:A=stages[-2]['a'];B=0
   else:A=pa^pb;B=pb
  b=v[i]>>12;a=(v[i]&4095)^F.mul(nodes[i]&4095,b);aa=F.mul(g,a)^A;bb=F.mul(g,b)^B
  if not aa and not bb:
   if all(bmul(v[j],g)^A^bmul(nodes[j],B)==0 for j in range(i,n)):break
   exc='noncompletion_zero_value';break
  prior_zero=bool(i and not stages[-1]['b'])
  kind='gap_second' if prior_zero else 'gap_first' if not bb else 'normal'
  for j in range(i+1,n):
   den=v[j]^a^bmul(nodes[j],b)
   if not den:exc='extension_pole';break
   v[j]=bmul(einv(den),kap[i]^kap[j])
  if exc:break
  previous_a=stages[-1]['a'] if i else 0
  stages.append(dict(kind=kind,a=aa,b=bb,kappa=kap[i],previous_a=previous_a))
  if kind=='gap_second' and not bb:exc='overlapping_scalar_chart';break
  factor=1 if not bb or prior_zero else F.mul(bb,bb);g=F.mul(factor,F.inv(g))
 if stages and stages[-1]['kind']=='gap_first':exc='unpaired_scalar_chart'
 return stages,exc


def synthesize(stages):
 N=[];B=[1]
 for s in reversed(stages):
  a,b,k=s['a'],s['b'],s['kappa'];q=[k,1,1]
  if s['kind']=='normal':
   h=F.inv(b);r=F.mul(a,h);c=F.mul(norm(a,b,k),h)
   N,B=F.add(F.product([r,1],N),F.scale(B,c)),F.add(F.scale(N,h),F.product([r^1,1],B))
  elif s['kind']=='gap_first':N,B=F.add(F.scale(N,a),F.product(q,B)),N
  else:
   olda=s['previous_a'];v=[a,b]
   N,B=F.add(F.product(v,N),F.product(F.add(q,F.scale(v,olda)),B)),F.add(N,F.scale(B,olda))
 scale=F.inv(B[-1]);return F.scale(N,scale),F.scale(B,scale)


def check(N=70,cases=720):
 start=time.time();rng=random.Random(142058+N);n=N//2;nodes,kap,P,maps=node_maps(n);exceptions={};counts={};passed=0;forced=0
 for trial in range(cases):
  d=trial%(n+1);pos=rng.sample(range(1,4096),d);val=[rng.randrange(1,4096) for _ in pos]
  if trial%2 and d>=2:
   coefficients=[F.mul(F.mul(a,F.evaluate(P,a)),F.inv(F.mul(a,a)^a^K)) for a in pos]
   x=0
   for c,e in zip(coefficients[:-1],val[:-1]):x^=F.mul(c,e)
   val[-1]=F.mul(x,F.inv(coefficients[-1]))
   if not val[-1]:continue
   forced+=1
  stages,exc=extract(N,F.syndrome(pos,val,N))
  if exc:exceptions[exc]=exceptions.get(exc,0)+1;continue
  nn,bb=synthesize(stages);want=[1]
  for a in pos:want=F.product(want,[a,1])
  assert bb==want,dict(trial=trial,degree=d,stages=stages,den=bb,want=want)
  for a,e in zip(pos,val):
   residue=F.mul(F.evaluate(nn,a),F.inv(F.evaluate(F.derivative(bb),a)));got=F.mul(residue,F.inv(F.mul(a,F.evaluate(P,a))));assert got==e,(trial,a,e,got)
  for s in stages:counts[s['kind']]=counts.get(s['kind'],0)+1
  passed+=1
 r=dict(status='MODEL_PASS',N=N,cases=cases,passed=passed,forced_scalar_pivots=forced,exceptions=exceptions,stage_kinds=counts,seconds=time.time()-start,scope='fixed-node raw extraction and explicit scalar/gap matrix synthesis, true locations and all magnitudes; completed tails stopped classically, not a quantum resource claim');Path(f'results/schur_gap_n{N}_model.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
if __name__=='__main__':
 for N,cases in [(6,216),(14,216),(70,720)]:check(N,cases)
