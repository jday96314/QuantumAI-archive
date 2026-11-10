"""Coherent field scaling through a routed zero word.

The inherited selected scale builds and reversibly clears two 12-bit factors
equal to one off-branch (48 Toffolis). Route the selected target through a
zero word instead: two controlled swaps cost 24 Toffolis and use 12 fewer
qubits. The ordinary field scale and its measurement eraser stay unchanged.
Aliasing factors retain the original implementation. Prior files are read-only.
"""
from contextlib import contextmanager
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
GAP=(HERE/'../../5/dqi-gap').resolve()
if str(GAP) not in sys.path:sys.path.append(str(GAP))
import quantum_control as original
prior_selected=original.selected

@contextmanager
def selected(m,E,flag):
 with prior_selected(m,E,flag) as (prior,ce):
  class Module:
   def __getattr__(self,name):return getattr(prior,name)
   def mon_slot(self,ignored,word,scale,eraser,product):
    if set(word)&(set(scale)|set(eraser)):
     return prior.mon_slot(ignored,word,scale,eraser,product)
    with m.opened(E,'routed_selected_scale',12) as routed:
     E.swaps(flag,list(zip(word,routed)))
     m.mon_slot(E,routed,scale,eraser,product)
     E.swaps(flag,list(zip(word,routed)))
  yield Module(),ce

def install():
 for module in list(sys.modules.values()):
  namespace=getattr(module,'__dict__',None)
  if namespace is not None and namespace.get('selected') is prior_selected:namespace['selected']=selected
