"""Build GF16 inverse-node parities directly in existing input wires."""
from contextlib import contextmanager
import inspect
import kickmix as km
import cubic_mask_phase_joint as joint
from cubic_trace_field import inverse_network


@contextmanager
def inverse4(b,source):
 net=inverse_network();live={}
 def apply(a,c,target,back=False):
  left={source[j] if j<4 else live[j] for j in a.terms};right={source[j] if j<4 else live[j] for j in c.terms};ca=a.const;cb=c.const
  if not left or not right:
   s=c if not left and ca else a if not right and cb else None
   if s is not None:
    if s.const:b.x(target)
    for j in s.terms:b.cx(source[j] if j<4 else live[j],target)
   return
  if left==right:
   if ca!=cb:return
   if ca:b.x(target)
   for q in left:b.cx(q,target)
   return
  p=min(left,key=lambda q:q.id);first=sorted(left-{p},key=lambda q:q.id)
  for q in first:b.cx(q,p)
  if p in right:right.symmetric_difference_update(left-{p})
  q=min(right-{p},key=lambda q:q.id);second=sorted(right-{q},key=lambda q:q.id)
  for r in second:b.cx(r,q)
  if ca:b.x(p)
  if cb:b.x(q)
  (b.del_and if back else b.init_and)(p,q,target)
  if cb:b.x(q)
  if ca:b.x(p)
  for r in reversed(second):b.cx(r,q)
  for r in reversed(first):b.cx(r,p)
 for j,(a,c) in enumerate(net.nodes):
  q=b.alloc_qubits(1)[0];apply(a,c,q);live[4+j]=q
 yield [(s.const,[source[j] if j<4 else live[j] for j in s.terms]) for s in net.outputs]
 for j in reversed(range(len(net.nodes))):
  a,c=net.nodes[j];q=live[4+j];apply(a,c,q,True);del live[4+j];b.free(km.array([q]))


def phase_pair(b,word,signals):
 for q,(const,terms) in zip(word,signals):
  if const:b.z(q)
  for r in terms:b.cz(q,r)
ns=dict(joint.ns);ns.update(vars(joint));ns.update(inverse4=inverse4,phase_pair=phase_pair);exec(joint.s.replace(joint.old,joint.new).replace('b.cz(Q,iv)','phase_pair(b,Q,iv)'),ns);exec(inspect.getsource(joint.phase_inverse),ns);phase_inverse=ns['phase_inverse'];exec(inspect.getsource(joint.native_phase_inverse),ns);native_phase_inverse=ns['native_phase_inverse']
if __name__=='__main__':
 exec(inspect.getsource(joint.test).replace('cubic_mask_phase_joint_native','cubic_mask_phase_lean_native'),ns);ns['test']()
