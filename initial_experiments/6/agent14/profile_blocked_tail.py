"""Exact resources for proposed blocked CF tail, excluding head/matrix cost.

This deliberately gives a lower bound on the candidate implementation, not a
complete circuit estimate; scalar pivots are arbitrary nonzero field words.
"""
import json
from pathlib import Path
from import_published import load
from low_workspace_inverse import bind
from triangular_dp import DPArithmetic
from profile_hybrid import DiscardTape
from cfraction_decoder import scale

def profile(n,b):
 m=load();m.inv_word=bind(m);B=m.Builder();B._items=DiscardTape();W=[B.qubits(12) for _ in range(70)]
 # CF history and tail occupy the original input bank. Local D,L,U,K
 # words are additional space in this optimistic isolated stage.
 t=W[70-n:];d=b//2
 E=m.Em(B,*B.qubits(3),m.Alloc(B));math=DPArithmetic(m,E)
 DL=[B.qubits(12) for _ in range(d+1)];L=[DL[0]]+[B.qubits(12) for _ in range(d-1)];U=[B.qubits(12) for _ in range(d)];K=B.qubits(12)
 math.triangular(DL,1,t)
 for s,w in zip(U,t):
  for q,r in zip(s,w):E.cx(q,r)
 scale(m,E,t,K,True)
 math.triangular(L,1,t)
 math.inverse([DL[0],*t])
 math.triangular(L,1,t,True)
 math.triangular(DL,1,t)
 return dict(n=n,block=b,Q=B.peak_qubits,T=B._cur_stats['toffoli'])

if __name__=='__main__':
 rows=[]
 for b in [2,4,6,8,10]:
  rr=[]
  for n in range(70-b,0,-b):
   r=profile(n,b);rr.append(r);print(json.dumps(r),flush=True)
  rows.append(dict(block=b,tail_forward_T=sum(r['T'] for r in rr),Q=max(r['Q'] for r in rr),rows=rr))
  Path('results/blocked_cf_tail_profiles.json').write_text(json.dumps(dict(scope='lower bound excludes head extraction, matrix formation, normalization product, synthesis, reversal and scan',families=rows),indent=2)+'\n')
